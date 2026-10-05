"""Build and register a private Pathway Tools image from a local installer."""
import argparse
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import uuid
from urllib.request import build_opener, HTTPRedirectHandler

from metapathways import nextflow


# Changing the recipe produces a distinct image, even for the same installer.
RECIPE = r'''Bootstrap: docker
From: ubuntu:22.04

%files
    installer /opt/pt-installer
    official-patches /opt/mp-official-patches

%post
    set -eu
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends ca-certificates xterm openssl libxml2 xvfb xauth libxm4 libssl-dev procps bzip2 ncbi-blast+
    mkdir -p /data /opt/bin
    chmod 700 /opt/pt-installer
    unset DISPLAY
    printf '/opt/pathway-tools\n/data\n\nn\nY\nn\n\n' | /opt/pt-installer
    test -x /opt/pathway-tools/pathway-tools
    rm /opt/pt-installer
    # Install only the unmodified, release-specific vendor files (SRI FAQ 7.4).
    pt_version=$(cat /opt/mp-official-patches/version)
    patch_dir=/opt/pathway-tools/aic-export/pathway-tools/ptools/$pt_version/patches
    test -d "$patch_dir"
    set -- "$patch_dir"/bin-*
    test "$#" -eq 1 && test -d "$1"
    for patch in /opt/mp-official-patches/files/*; do
        case "$patch" in
            *.fasl) cp "$patch" "$1/" ;;
            *) cp "$patch" "$patch_dir/" ;;
        esac
    done
    # Load official patches while the filesystem is still writable.
    xvfb-run -a /opt/pathway-tools/pathway-tools -no-patch-download -lisp -eval '(progn (format t "~%MP-PT-READY~%") (exit))' > /opt/mp-pt-patch-startup.log 2>&1 || { cat /opt/mp-pt-patch-startup.log; exit 1; }
    cat /opt/mp-pt-patch-startup.log
    grep -qx 'MP-PT-READY' /opt/mp-pt-patch-startup.log
    printf '[ncbi]\nData=/usr/share/ncbi/data\n' > /opt/mp-ncbirc
    test -d /usr/share/ncbi/data
    blastp -version
    makeblastdb -version
    dpkg-query -W ncbi-blast+ > /opt/mp-blast-version.txt
    # Keep a pristine template; each invocation mounts its own writable /data.
    cp -a /data/ptools-local /opt/ptools-local-template
    rm -rf /var/lib/apt/lists/*

%environment
    export PATH=/opt/pathway-tools:$PATH

%runscript
    exec /opt/pathway-tools/pathway-tools "$@"

%labels
    org.metapathways.purpose PathwayTools
'''


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError('Official patch URL redirected; refusing an unverified patch source: ' + newurl)


class _PatchLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.names = set()

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            href = dict(attrs).get('href', '')
            # No arbitrary URLs, parent paths, nested directories or query strings.
            if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.(?:fasl|lisp|tar\.gz)', href):
                self.names.add(href)


def installer_version(installer, explicit=None):
    match = re.search(r'pathway-tools-(\d+\.\d+)', Path(installer).name)
    version = explicit or (match.group(1) if match else None)
    if not version or not re.fullmatch(r'\d+\.\d+', version):
        raise ValueError('Cannot determine Pathway Tools release; use --ptools_version (for example 29.5) for a renamed installer')
    if explicit and match and explicit != match.group(1):
        raise ValueError('--ptools_version disagrees with the installer filename')
    return version


def download_patches(version, destination):
    """Snapshot SRI's official release feed; never silently use partial downloads."""
    if not re.fullmatch(r'\d+\.\d+', version):
        raise ValueError('Invalid Pathway Tools release')
    url = f'https://bioinformatics.ai.sri.com/ptools/{version}/Linux-64/patches/'
    destination = Path(destination)
    files = destination / 'files'
    files.mkdir(parents=True)
    opener = build_opener(_NoRedirect())
    try:
        with opener.open(url, timeout=60) as response:
            listing = response.read()
        links = _PatchLinks()
        links.feed(listing.decode('utf-8'))
        if not links.names:
            raise ValueError('Vendor listing contained no recognized patch files')
        manifest = dict(source=url, version=version, fetched_at=datetime.now(timezone.utc).isoformat(),
                        listing_sha256=hashlib.sha256(listing).hexdigest(), files=[])
        for name in sorted(links.names):
            print('Downloading official Pathway Tools patch: ' + name, flush=True)
            target = files / name
            with opener.open(url + name, timeout=60) as response, target.open('wb') as stream:
                shutil.copyfileobj(response, stream)
            if target.stat().st_size == 0:
                raise ValueError('Empty patch: ' + name)
            manifest['files'].append(dict(name=name, url=url + name, sha256=digest(target)))
        manifest['snapshot_sha256'] = hashlib.sha256(json.dumps(manifest['files'], sort_keys=True).encode()).hexdigest()
        (destination / 'version').write_text(version + '\n')
        (destination / 'index.html').write_bytes(listing)
        save_json(destination / 'manifest.json', manifest)
        return manifest
    except Exception as exc:
        raise RuntimeError(f'Official Pathway Tools patches could not be obtained from {url}; build stopped: {exc}') from exc


def registry_path():
    return Path(os.environ.get('XDG_CONFIG_HOME', Path.home() / '.config')) / 'metapathways' / 'ptools.json'


def registered_image():
    override = os.environ.get('METAPATHWAYS_PTOOLS_IMAGE')
    if override:
        image = Path(override).expanduser().resolve()
    elif registry_path().exists():
        image = Path(json.loads(registry_path().read_text())['image'])
    else:
        return None
    if not image.is_file():
        raise ValueError(f'Registered Pathway Tools image is missing: {image}; run metapathways build_pt again')
    return str(image)


def save_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        temp = Path(f.name)
        json.dump(data, f, indent=2)
        f.write('\n')
    try:
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def exec_command(image, state, command):
    """Use a task-private data directory, PID namespace, home and /tmp."""
    state = Path(state).resolve()
    state.mkdir(parents=True, exist_ok=True)
    if any(c in str(state) for c in ':,\n'):
        raise ValueError('Apptainer state path cannot contain a colon, comma or newline')
    executable = shutil.which('apptainer')
    if not executable:
        raise RuntimeError('Apptainer is required and must be on PATH')
    return [executable, 'exec', '--containall', '--cleanenv',
            '--home', f'{state}:/data', '--pwd', '/data', str(image), *command]


def validate(image, directory):
    script = ('cp -a /opt/ptools-local-template /data/ptools-local; '
              'cp /opt/mp-ncbirc /data/.ncbirc; '
              'blastp -version; makeblastdb -version; '
              'printf ">mp_validation\\nMKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE\\n" > /data/check.faa; '
              'makeblastdb -in /data/check.faa -dbtype prot -out /data/check-db; '
              'blastp -query /data/check.faa -db /data/check-db -outfmt 6 -out /data/check.tsv; '
              'test -s /data/check.tsv; '
              'exec xvfb-run -a /opt/pathway-tools/pathway-tools '
              '-no-patch-download -lisp -eval \'(progn (format t "~%MP-PT-READY~%") (exit))\'')
    result = subprocess.run(exec_command(image, directory, ['sh', '-ec', script]),
                            check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=300)
    if result.returncode:
        raise RuntimeError(f'Pathway Tools startup failed (exit {result.returncode}): ' + result.stdout[-4000:])
    if 'MP-PT-READY' not in [line.strip() for line in result.stdout.splitlines()]:
        raise RuntimeError('Pathway Tools startup did not emit its validation marker: ' + result.stdout[-4000:])
    if re.search(r'Error determining path to BLAST|blastall or blastp could not be located', result.stdout, re.I):
        raise RuntimeError('Pathway Tools could not locate the installed BLAST executables: ' + result.stdout[-4000:])
    return result.stdout


def run_pgdb(image, inputs, outputs, tag, taxprune=True, sample_output=None, transport_inference=True):
    """Run one PGDB with no shared host Pathway Tools state."""
    if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_-]*', tag):
        raise ValueError('PGDB tag must start with a letter or underscore and contain only letters, digits, underscores or hyphens')
    output = Path(outputs).resolve()
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.pt-run-', dir=output.parent) as work:
        state = Path(work)
        shutil.copytree(inputs, state / 'input')
        (state / 'output').mkdir()
        script = '''set -eu
cp -a /opt/ptools-local-template /data/ptools-local
if test -f /opt/mp-ncbirc; then cp /opt/mp-ncbirc /data/.ncbirc; fi
mkdir -p /data/blastdb /tmp/.X11-unix
# --containall isolates /tmp, but not Linux abstract sockets. Use only the
# private filesystem X socket so concurrent containers cannot collide.
printf 'build\n' > /data/stage.txt
xvfb-run -a -e /data/build-xvfb.log -s '-screen 0 1280x1024x24 -nolisten local' /opt/pathway-tools/pathway-tools -patho /data/input "$@" -no-web-cel-overview -no-patch-download -no-cel-overview -disable-metadata-saving -nologfile
'''
        lisp = f"(progn (with-organism (:org-id '{tag}) (dump-frames-to-attribute-value-files (org-data-dir)))(exit))"
        script += "printf 'export\\n' > /data/stage.txt\n"
        script += shlex.join(['xvfb-run', '-a', '-e', '/data/export-xvfb.log',
                              '-s', '-screen 0 1280x1024x24 -nolisten local', '/opt/pathway-tools/pathway-tools',
                              '-no-patch-download', '-no-cel-overview', '-nologfile', '-eval', lisp]) + '\n'
        script += "printf 'archive\\n' > /data/stage.txt\n"
        script += shlex.join(['tar', '-cjf', f'/data/output/{tag}cyc.tar.bz2',
                              '-C', f'/data/ptools-local/pgdbs/user/{tag.lower()}cyc', '.']) + '\n'
        command = ['sh', '-ec', script, 'run-pgdb']
        if transport_inference:
            command.append('-tip')
        if not taxprune:
            command.append('-no-taxonomic-pruning')
        invocation = exec_command(image, state, command)
        diagnostics = output / 'diagnostics' / uuid.uuid4().hex
        status = dict(tag=tag, image=str(image), command=invocation, status='FAILED')
        try:
            if sample_output is not None:
                from metapathways.pt_sequences import attach_sequences
                (state / 'stage.txt').write_text('input preparation\n')
                attach_sequences(state / 'input', sample_output)
            from metapathways.pt_reactions import filter_reactions
            filter_reactions(state / 'input', image=image, sample_output=sample_output)
            subprocess.run(invocation, check=True)
            archive = state / 'output' / f'{tag}cyc.tar.bz2'
            if not archive.is_file() or not archive.stat().st_size:
                raise RuntimeError(f'Pathway Tools produced no PGDB archive for {tag}')
            archive.replace(output / archive.name)
            status['status'] = 'SUCCESS'
        except BaseException as exc:
            status['error'] = str(exc)
            status['exit_code'] = getattr(exc, 'returncode', None)
            # Pathologic redirects its own errors away from the parent's stderr.
            private_log = state / 'input/pathologic.log'
            if private_log.is_file():
                from collections import deque
                with private_log.open(errors='replace') as stream:
                    tail = ''.join(deque(stream, maxlen=80))
                print(f'Pathway Tools internal log tail ({tag}):\n{tail}', flush=True)
            raise
        finally:
            # Preserve internal diagnostics before TemporaryDirectory removes state,
            # on success, failure, or an interrupted subprocess.
            diagnostics.mkdir(parents=True, exist_ok=True)
            stage = state / 'stage.txt'
            status['stage'] = stage.read_text().strip() if stage.is_file() else 'container startup'
            for source in state.rglob('*'):
                if source.is_file() and (source.suffix.lower() in ('.log', '.err', '.out')
                                         or source.name in ('stage.txt', 'sequence-input.json', 'ptools-reaction-filter.json', 'organism-params.dat')
                                         or 'reports' in source.relative_to(state).parts):
                    destination = diagnostics / source.relative_to(state)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, destination)
            if status['status'] != 'SUCCESS':
                # Retain the last on-disk PGDB for separate recovery attempts.
                # Moving avoids copying a potentially large database before cleanup.
                user_pgdbs = state / 'ptools-local/pgdbs/user'
                if user_pgdbs.is_dir():
                    retained = diagnostics / 'failed-pgdbs'
                    shutil.move(str(user_pgdbs), str(retained))
                    status['retained_pgdbs'] = str(retained)
                shutil.copytree(state / 'input', diagnostics / 'input', dirs_exist_ok=True)
            save_json(diagnostics / 'execution.json', status)
            print(f'Pathway Tools {status["status"]} during {status["stage"]}; diagnostics: {diagnostics}', flush=True)


def build(installer, image, installer_sha256, threads, version=None):
    """Worker entry point: never publish or register a failed build."""
    installer, image = Path(installer).resolve(), Path(image).resolve()
    if digest(installer) != installer_sha256:
        raise ValueError('Installer changed after the build was planned')
    executable = shutil.which('apptainer')
    if not executable:
        raise RuntimeError('Apptainer is required and must be on PATH')
    image.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.pt-build-', dir=image.parent) as work:
        work = Path(work)
        # A fixed relative name avoids definition-file quoting of arbitrary paths.
        (work / 'installer').symlink_to(installer)
        patches = download_patches(installer_version(installer, version), work / 'official-patches')
        definition = work / 'pathway-tools.def'
        definition.write_text(RECIPE)
        partial = work / 'pathway-tools.sif'
        subprocess.run([executable, 'build', '--fakeroot', '--mksquashfs-args',
                        f'-processors {threads}', str(partial), str(definition)],
                       cwd=work, check=True)
        validation_log = validate(partial, work / 'validation')
        metadata = dict(installer_sha256=installer_sha256, image=str(image),
                        image_sha256=digest(partial), recipe_sha256=hashlib.sha256(RECIPE.encode()).hexdigest(),
                        official_patches=patches, validation_log=validation_log)
        partial.replace(image)
        save_json(str(image) + '.json', metadata)
        image.with_suffix('.def').write_text(RECIPE)


def parser():
    p = argparse.ArgumentParser(prog='metapathways build_pt',
                                description='Build a Pathway Tools SIF from a local Linux installer using Nextflow and Apptainer.')
    p.add_argument('-i', '--installer', required=True, help='local Pathway Tools Linux x86-64 installer')
    p.add_argument('--ptools_version', help='release number if the installer was renamed; otherwise inferred from its filename')
    p.add_argument('-d', '--refdb_dir', help='also export and prepare the licensed MetaCyc reference in this MPDB after building the SIF')
    p.add_argument('--skip_pt_screen', action='store_true', help='Skip default reaction compatibility screening when building MetaCyc with -d')
    p.add_argument('-a', '--aligner', choices=('fast', 'blast'), default='fast', help='MetaCyc reference index format with -d [fast]')
    default = Path(os.environ.get('XDG_DATA_HOME', Path.home() / '.local/share')) / 'metapathways/containers'
    p.add_argument('-o', '--output_dir', default=str(default), help='container directory [~/.local/share/metapathways/containers]')
    p.add_argument('-t', '--threads', type=nextflow.positive, default=2, help='CPUs for image compression [2]; Pathway Tools uses one CPU')
    p.add_argument('--dryrun', action='store_true', help='write and show the build plan without building or registering')
    nextflow.add_resources(p, task_memory='4 GB')
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    installer = Path(args.installer).expanduser().resolve()
    if not installer.is_file():
        p.error(f'Installer does not exist: {installer}')
    if platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'AMD64'):
        p.error('Building this Pathway Tools image requires Linux x86-64')
    if not args.dryrun:
        for executable in ('apptainer', 'nextflow'):
            if not shutil.which(executable):
                p.error(f'{executable} must be installed and on PATH')
    checksum = digest(installer)
    recipe_hash = hashlib.sha256(RECIPE.encode()).hexdigest()
    try:
        version = installer_version(installer, args.ptools_version)
    except ValueError as exc:
        p.error(str(exc))
    output = Path(args.output_dir).expanduser().resolve()
    # Each explicit build refreshes the vendor feed and creates a separate SIF.
    # Never overwrite an image that existing analysis receipts/running jobs use.
    build_id = uuid.uuid4().hex
    image = output / f'pathway-tools-{version}-{checksum[:12]}-{recipe_hash[:8]}-{build_id[:12]}.sif'
    command = shlex.join([sys.executable, '-m', 'metapathways.pt_container', '--worker',
                          str(installer), str(image), checksum, str(args.threads), version])
    t = nextflow.task('build-pt-' + checksum + recipe_hash + build_id, 'Build Pathway Tools ' + version,
                     [command], [str(installer)], [str(image), str(image) + '.json'],
                     cpus=args.threads, memory=args.memory, adopt_existing=False)
    tasks = [t]
    if args.refdb_dir:
        from metapathways.nf_databases import metacyc_task, screen_task
        tasks.append(metacyc_task(args.refdb_dir, image, args.aligner, args.memory, [t['id']]))
        if not args.skip_pt_screen:
            tasks.append(screen_task(args.refdb_dir, image, args.memory, resources=args))
    nextflow.launch(tasks, output, args, 'build_pt', dryrun=args.dryrun)
    if not args.dryrun:
        metadata = json.loads(Path(str(image) + '.json').read_text())
        if metadata['image_sha256'] != digest(image):
            raise RuntimeError('Built image checksum does not match its metadata')
        save_json(registry_path(), metadata)
        print(f'Pathway Tools image: {image}\nRegistered in: {registry_path()}')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        build(*sys.argv[2:5], int(sys.argv[5]), sys.argv[6])
    else:
        main()
