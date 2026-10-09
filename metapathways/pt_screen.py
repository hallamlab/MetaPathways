"""Maintainer compatibility screen for explicit MetaCyc reaction assignments."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
import time

from metapathways.pt_container import exec_command, registered_image
from metapathways.nf_worker import atomic_json

SCREEN_VERSION = 1


def digest(file):
    h = hashlib.sha256()
    with Path(file).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def reactions_from_db(directory):
    file = Path(directory)/'functional_categories/MetaCyc-monomer-rxn-pairs.tsv'
    with file.open() as stream:
        rows = csv.DictReader(stream, delimiter='\t')
        if not {'MC', 'RXN'} <= set(rows.fieldnames or []):
            raise ValueError(f'Missing MC/RXN columns: {file}')
        reactions = sorted({row['RXN'].strip() for row in rows if row['RXN'].strip()})
    if not reactions:
        raise ValueError(f'No reaction assignments found: {file}')
    if any(not re.fullmatch(r'[A-Za-z0-9_.:+-]+', r) for r in reactions):
        raise ValueError('Unsafe reaction identifier in mapping table')
    return reactions, file


def write_inputs(directory, reactions, tag):
    directory.mkdir(parents=True)
    count = max(1, len(reactions))
    gene = 'ATG' + 'GCT'*78 + 'TAA'
    sequence = gene * count
    (directory/'contig.fasta').write_text('>contig\n'+sequence+'\n')
    records = []
    for i in range(count):
        records.append(f'ID\tG{i+1}\nNAME\tG{i+1}\nSTARTBASE\t{i*len(gene)+1}\n'
                       f'ENDBASE\t{(i+1)*len(gene)}\nFUNCTION\tUncharacterized protein\n'
                       'PRODUCT-TYPE\tP\n' +
                       (f'METACYC\t{reactions[i]}\n' if reactions else '') + '//\n')
    (directory/'contig.pf').write_text(''.join(records))
    (directory/'genetic-elements.dat').write_text('ID\tcontig\nNAME\tcontig\nTYPE\t:CONTIG\n'
        'CODON-TABLE\t11\nANNOT-FILE\tcontig.pf\nSEQ-FILE\tcontig.fasta\n//\n')
    (directory/'organism-params.dat').write_text(f'ID\t{tag}\nSTORAGE\tFILE\nNAME\t{tag}\n'
        'STRAIN\t1\nRANK\t|species|\nNCBI-TAXON-ID\t131567\n')


def execute_attempt(image, directory, reactions, timeout, scratch=None):
    """Use unfiltered synthetic inputs; preserve logs and a real success marker."""
    tag = 'MPscreen'
    write_inputs(directory/'input', reactions, tag)
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='mp-pt-screen-', dir=scratch) as temporary:
        state = Path(temporary)
        # Bind only this attempt's inputs into its private writable home.
        import shutil
        shutil.copytree(directory/'input', state/'input')
        (state/'screen.lisp').write_text("""(progn
 (handler-case
  (progn
   (batch-pathologic "1.0" "/data/input/"
      :download-publications? nil :do-overview? nil :web-cel-ov? nil
      :taxonomic-pruning? t :tip? t :suppress-metadata-saving? t
      :standard-streams? t :debug? t :trap-errors? nil)
   (format t "~%MP-SCREEN-SUCCESS~%") (finish-output) (exit))
  (error (e) (format t "~%MP-SCREEN-ERROR ~A~%" e)
             (finish-output) (excl:exit 1))))
""")
        script = """set -eu
cp -a /opt/ptools-local-template /data/ptools-local
if test -f /opt/mp-ncbirc; then cp /opt/mp-ncbirc /data/.ncbirc; fi
exec xvfb-run -a -s '-screen 0 1280x1024x24 -nolisten local' /opt/pathway-tools/pathway-tools -no-patch-download -no-cel-overview -nologfile -lisp -load /data/screen.lisp
"""
        with (directory/'console.log').open('w') as log:
            process = subprocess.Popen(exec_command(image, state, ['sh', '-ec', script]),
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            timed_out = False
            try:
                process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
        for source in state.rglob('*.log'):
            destination = directory/'diagnostics'/source.relative_to(state)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    text = (directory/'console.log').read_text(errors='replace')
    if timed_out or process.returncode < 0 or process.returncode in (137, 143):
        status = 'INCONCLUSIVE'
    elif process.returncode == 0 and 'MP-SCREEN-SUCCESS' in text.splitlines():
        status = 'PASS'
    elif 'MP-SCREEN-ERROR' in text:
        status = 'FAIL'
    else:
        status = 'INCONCLUSIVE'
    return dict(status=status, exit_code=process.returncode, timed_out=timed_out,
                duration_seconds=round(time.monotonic()-start, 3))


class Screen:
    def __init__(self, args, runner=execute_attempt):
        self.args, self.runner = args, runner
        self.output = Path(args.output_dir).expanduser().resolve()
        if getattr(args, 'database_screen', False):
            mapping = Path(args.refdb_dir)/'functional_categories/MetaCyc-monomer-rxn-pairs.tsv'
            self.output /= digest(mapping)
        self.output.mkdir(parents=True, exist_ok=True)
        self.image = Path(args.image).expanduser().resolve()
        self.reactions, mapping = reactions_from_db(args.refdb_dir)
        if args.reactions:
            missing = set(args.reactions)-set(self.reactions)
            if missing:
                raise ValueError('Reactions absent from MPDB: '+', '.join(sorted(missing)))
            self.reactions = sorted(set(args.reactions))
        self.identity = dict(screen_version=SCREEN_VERSION, image_sha256=digest(self.image),
            mapping_sha256=digest(mapping), reactions=self.reactions, batch_size=args.batch_size,
            confirm_runs=args.confirm_runs, timeout=args.timeout)
        manifest = self.output/'screen.json'
        if manifest.exists() and json.loads(manifest.read_text()) != self.identity:
            raise ValueError('Screen inputs/settings changed; use a new output directory')
        atomic_json(manifest, self.identity)
        self.candidates, self.interactions, self.inconclusive = {}, [], []

    def attempt(self, reactions, label):
        key = hashlib.sha256(json.dumps([reactions, label]).encode()).hexdigest()[:20]
        directory = self.output/'attempts'/key
        receipt = directory/'result.json'
        if receipt.exists():
            previous = json.loads(receipt.read_text())
            if previous['status'] != 'INCONCLUSIVE':
                return previous
        if directory.exists():
            # An interrupted attempt has no valid receipt; preserve its evidence.
            directory.rename(directory.with_name(key+'-interrupted-'+str(time.time_ns())))
        directory.mkdir(parents=True)
        print(f'Screen {label}: {len(reactions)} reactions; {directory}', flush=True)
        result = self.runner(self.image, directory, reactions, self.args.timeout, self.args.scratch_dir)
        result.update(reactions=reactions, label=label, directory=str(directory))
        atomic_json(receipt, result)
        print(f'Screen {label}: {result["status"]}', flush=True)
        return result

    def isolate(self, reactions):
        result = self.attempt(reactions, 'screen')
        if result['status'] == 'PASS':
            return []
        if result['status'] != 'FAIL':
            self.inconclusive.append(result)
            return []
        if len(reactions) > 1:
            half = len(reactions)//2
            failures = self.isolate(reactions[:half])+self.isolate(reactions[half:])
            if not failures:
                children = []
                for child in (reactions[:half], reactions[half:]):
                    key = hashlib.sha256(json.dumps([child, 'screen']).encode()).hexdigest()[:20]
                    children.append(json.loads((self.output/'attempts'/key/'result.json').read_text()))
                if all(child['status'] == 'PASS' for child in children):
                    self.interactions.append(result)
                else:
                    self.inconclusive.append(result)
            return failures
        checks = [result]
        for i in range(1, self.args.confirm_runs):
            checks.append(self.attempt(reactions, f'confirmation-{i}'))
        control = self.attempt([], 'control-'+reactions[0])
        if all(r['status'] == 'FAIL' for r in checks) and control['status'] == 'PASS':
            self.candidates[reactions[0]] = dict(reason='Repeated isolated PTools build failure; no-reaction control passed',
                evidence=[r['directory'] for r in checks], control=control['directory'],
                image_sha256=self.identity['image_sha256'], review_required=True)
            return reactions
        self.inconclusive.extend(checks+[control])
        return []

    def run(self):
        if self.attempt([], 'baseline')['status'] != 'PASS':
            raise RuntimeError('Baseline PGDB failed; inspect its logs before screening reactions')
        batches = [self.reactions[i:i+self.args.batch_size]
                   for i in range(0, len(self.reactions), self.args.batch_size)]
        # Every attempt has its own home and receipt. Results are consolidated after jobs finish.
        with ThreadPoolExecutor(max_workers=self.args.max_tasks) as pool:
            list(pool.map(self.isolate, batches))
        atomic_json(self.output/'blacklist-candidates.json', self.candidates)
        atomic_json(self.output/'summary.json', dict(reactions=len(self.reactions),
            initial_batches=len(batches), candidates=len(self.candidates),
            interaction_failures=self.interactions, inconclusive=self.inconclusive,
            image_sha256=self.identity['image_sha256'],
            note='Explicit-ID synthetic screen; passing does not certify all annotations or reaction combinations.'))
        if getattr(self.args, 'publish', False):
            if self.inconclusive or self.interactions:
                raise RuntimeError('Screen is unresolved; compatibility list was not published. Inspect summary.json and resume.')
            atomic_json(Path(self.args.refdb_dir)/'functional_categories/ptools_reaction_compatibility.json',
                        dict(image_sha256=self.identity['image_sha256'],
                             mapping_sha256=self.identity['mapping_sha256'], reactions=self.candidates,
                             screen_version=SCREEN_VERSION))
        print(f'Screen finished: {len(self.candidates)} blacklist candidates; '
              f'{len(self.interactions)} batch interaction failures; '
              f'{len(self.inconclusive)} inconclusive results. Review {self.output}', flush=True)


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return value


def main(argv=None):
    p = argparse.ArgumentParser(description='Screen explicit MPDB reaction assignments against a licensed PTools SIF')
    p.add_argument('-d', '--refdb_dir', required=True)
    p.add_argument('-o', '--output_dir', required=True)
    p.add_argument('--image', default=registered_image(), help='PTools SIF [registered build_pt image]')
    p.add_argument('--database_screen', action='store_true', help=argparse.SUPPRESS)
    p.add_argument('--publish', action='store_true', help='Save a completed full screen compatibility list in the MPDB')
    p.add_argument('--reactions', nargs='+', help='Optional reaction IDs for a targeted screen')
    p.add_argument('--batch_size', type=positive, default=100)
    p.add_argument('--max_tasks', type=positive, default=1, help='Concurrent isolated containers [1]')
    p.add_argument('--confirm_runs', type=positive, default=2)
    p.add_argument('--timeout', type=positive, default=1800, help='Seconds per build [1800]; timeout is inconclusive')
    p.add_argument('--scratch_dir', help='Local temporary storage for private PGDB builds')
    args = p.parse_args(argv)
    if not args.image or not Path(args.image).expanduser().is_file():
        p.error('Provide a valid --image or register one with build_pt')
    if args.publish and args.reactions:
        p.error('--publish requires the complete MPDB reaction set')
    if args.confirm_runs < 2:
        p.error('--confirm_runs must be at least 2')
    try:
        Screen(args).run()
    except (OSError, ValueError, RuntimeError) as exc:
        p.exit(1, f'screen_pt: {exc}\n')


if __name__ == '__main__':
    main()
