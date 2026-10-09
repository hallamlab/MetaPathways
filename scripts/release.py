#!/usr/bin/env python3
"""Prepare, build, validate, and publish MetaPathways releases.

Run with --help. 'publish', 'github-release', and 'upload-conda' write to remote services.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "hallamlab/MetaPathways"
VERSION_RE = re.compile(r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)(rc[1-9][0-9]*)?")
STAGES = {
    "PREPROCESS_INPUT", "ORF_PREDICTION", "ORF_TO_AMINO", "FILTER_AMINOS",
    "FUNC_SEARCH:swissprot_test", "COMPUTE_REFSCORES",
    "PARSE_FUNC_SEARCH:swissprot_test", "SCAN_rRNA:barrnap",
    "SCAN_rRNA:SILVA_SSU_test", "SCAN_rRNA:SILVA_LSU_test", "SCAN_tRNA",
    "ANNOTATE_ORFS", "CREATE_ANNOT_REPORTS", "GENBANK_FILE",
    "PATHOLOGIC_INPUT", "COMPUTE_TPM",
}
OUTPUTS = [
    "genbank/k12_test.annot.gff", "genbank/k12_test.gbk",
    "results/annotation_table/k12_test.ORF_annotation_table.txt",
    "results/annotation_table/k12_test.functional_and_taxonomic_table.txt",
    "results/rpkm/k12_test.contig_counts.tsv", "ptools/0.pf",
]


def run(*args, cwd=ROOT, capture=False, log=None, env=None):
    args = [str(a) for a in args]
    # Secrets must be supplied through environment variables, never arguments.
    print("+ " + " ".join(args), flush=True)
    if log:
        print(f"  log: {log}", flush=True)
        with Path(log).open("w") as stream:
            subprocess.run(args, cwd=cwd, check=True, stdout=stream,
                           stderr=subprocess.STDOUT, env=env)
        return ""
    result = subprocess.run(args, cwd=cwd, check=True, text=True,
                            stdout=subprocess.PIPE if capture else None, env=env)
    return result.stdout.strip() if capture else ""


def validate_runtime_dependencies(metadata):
    """Compare runtime metadata using tooling from the release environment."""
    from packaging.version import Version

    if metadata["pip_present"]:
        raise ValueError("pip must remain build-time only in the Conda package")
    if tuple(metadata["python"]) < (3, 11):
        raise ValueError("Python >=3.11 required")
    for package, minimum in (("urllib3", "2.8.0"), ("setuptools", "83.0.0")):
        if Version(metadata["versions"][package]) < Version(minimum):
            raise ValueError(f"{package} >={minimum} required")


def version(root=ROOT):
    return version_from_text((root / "metapathways/_version.py").read_text())


def version_from_text(text):
    tree = ast.parse(text)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "__version__" for t in node.targets
        ):
            return valid_version(ast.literal_eval(node.value))
    raise ValueError("No __version__ found")


def valid_version(value):
    value = value.removeprefix("v")
    if not VERSION_RE.fullmatch(value):
        raise ValueError("Use X.Y.Z or X.Y.ZrcN, e.g. 3.5.0 or 3.5.1rc1; no .dev suffix.")
    return value


def recipe_build(text):
    match = re.search(r"^  number: ([0-9]+)$", text, re.M)
    if not match:
        raise ValueError("Recipe must declare an integer build number.")
    return int(match.group(1))


def release_tag(value, number):
    return f"v{value}" + (f"-build{number}" if number else "")


def tag_version(tag):
    match = re.fullmatch(r"v(.+?)(?:-build([1-9][0-9]*))?", tag)
    if not match:
        raise ValueError("Expected a release tag such as v3.5.1 or v3.5.1-build1.")
    return valid_version(match.group(1)), int(match.group(2) or 0)


def clean():
    if run("git", "status", "--porcelain", "--untracked-files=normal", capture=True):
        raise ValueError("Commit or stash your changes first; releases use committed source only.")


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(stream.read()).hexdigest()


def prepare(args):
    value = valid_version(args.version)
    recipe = ROOT / "conda_recipe/meta_template.yaml"
    match = re.search(r"^  number: ([0-9]+)$", recipe.read_text(), re.M)
    if not match:
        raise ValueError("Recipe must declare an integer build number.")
    build_number = args.build_number
    if build_number is None:
        build_number = int(match.group(1)) if value == version(ROOT) else 0
    if build_number < 0:
        raise ValueError("Build number must be nonnegative.")
    target = ROOT / "metapathways/_version.py"
    text = re.sub(r'^__version__ = .+$', f'__version__ = "{value}"',
                  target.read_text(), flags=re.M)
    text = re.sub(r'^__status__ = .+$',
                  '__status__ = "' + ("Release Candidate" if "rc" in value else "Release") + '"',
                  text, flags=re.M)
    target.write_text(text)
    recipe.write_text(re.sub(r"^  number: [0-9]+$", f"  number: {build_number}",
                             recipe.read_text(), flags=re.M))
    readme = ROOT / "README.md"
    readme.write_text(re.sub(r"https://img.shields.io/badge/Version-[^)]*",
                            f"https://img.shields.io/badge/Version-{value}-blue.svg",
                            readme.read_text()))
    citation = ROOT / "CITATION.cff"
    if citation.exists():
        text = re.sub(r"^version:.*\n?", "", citation.read_text(), flags=re.M)
        citation.write_text(text.rstrip() + f'\nversion: "{value}"\n')
    print(f"Prepared {value}, Conda build {build_number}. Review and commit before publishing.")


def validate_run(directory):
    sample = Path(directory) / "test/k12_test"
    text = (sample / "metapathways_steps_log.txt").read_text()
    successes = re.findall(r"^([^\t\n]+)\tSUCCESS\b", text, re.M)
    if set(successes) != STAGES or len(successes) != len(STAGES):
        raise ValueError(f"Integration test incomplete: expected {len(STAGES)} distinct successful stages; "
                         f"missing {sorted(STAGES - set(successes))}")
    for name in OUTPUTS:
        path = sample / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty output: {path}")
    for path in [sample / "errors_warnings_log.txt",
                 Path(directory) / "test/global_errors_warnings.txt"]:
        text = path.read_text()
        if re.search(r"(?im)^\s*(?:ERROR\b|FAILED\b|Traceback\b)", text):
            raise ValueError(f"Review errors in {path}")
    return {"successful_stages": sorted(STAGES), "nonempty_outputs": OUTPUTS}


def validate_citation(data, value):
    if not isinstance(data, dict) or data.get('version') != value:
        raise ValueError('CITATION.cff must match the release version.')
    authors = data.get('authors')
    if not isinstance(authors, list) or not authors:
        raise ValueError('CITATION.cff must contain the reviewed manuscript authors.')
    names = set()
    for author in authors:
        if not isinstance(author, dict) or not all(
            isinstance(author.get(key), str) and author[key].strip()
            for key in ('family-names', 'given-names', 'affiliation')
        ):
            raise ValueError('Each citation author needs a name and reviewed affiliation.')
        name = (author['family-names'].strip().casefold(), author['given-names'].strip().casefold())
        if name in names:
            raise ValueError('Duplicate citation author; review CITATION.cff.')
        names.add(name)


def build(args):
    clean()
    value = version()
    import yaml
    citation = ROOT / 'CITATION.cff'
    if not citation.is_file() or (ROOT / '.zenodo.json').exists():
        raise ValueError('Release requires CITATION.cff without an overriding .zenodo.json.')
    validate_citation(yaml.safe_load(citation.read_text()), value)
    number = recipe_build((ROOT / "conda_recipe/meta_template.yaml").read_text())
    tag = release_tag(value, number)
    if args.tag and args.tag != tag:
        raise ValueError(f"Tag {args.tag!r} does not match source version v{value}")
    commit = run("git", "rev-parse", "HEAD", capture=True)
    if args.tag and run("git", "rev-parse", f"{args.tag}^{{commit}}", capture=True) != commit:
        raise ValueError("Release tag must point to the checked-out commit.")
    if not args.source_only and (platform.system(), platform.machine()) != ("Linux", "x86_64"):
        raise ValueError("Bundled binaries require Linux x86-64.")
    output = Path(args.output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"{output} is not empty. Choose a fresh --output directory.")
    output.mkdir(parents=True, exist_ok=True)
    archive = output / f"metapathways-{value}-source.zip"
    run("git", "archive", "--format=zip", f"--prefix=metapathways-{value}/",
        "-o", archive, "HEAD")
    validation = {"scope": "source-only", "version": value, "commit": commit}
    with tempfile.TemporaryDirectory(prefix="metapathways-release-") as temp:
        work = Path(temp)
        with zipfile.ZipFile(archive) as source:
            source.extractall(work)
        snapshot = work / f"metapathways-{value}"
        # Restore executable bits: ZipFile.extractall does not preserve Git modes.
        with zipfile.ZipFile(archive) as source:
            for member in source.infolist():
                mode = member.external_attr >> 16
                if mode and not member.is_dir():
                    (work / member.filename).chmod(mode & 0o777)
        run(sys.executable, "-m", "build", "--sdist", "--no-isolation",
            "--outdir", output, cwd=snapshot, log=output / "source-build.log")
        sdist = output / f"metapathways-{value}.tar.gz"
        if not sdist.is_file():
            raise ValueError(f"Expected source distribution {sdist}")
        if not args.source_only:
            recipe = work / "recipe"
            run(sys.executable, snapshot / "conda_recipe/compile_recipe.py",
                "--sdist", sdist, "--output-dir", recipe, cwd=snapshot)
            channel = work / "channel"
            channel.mkdir()
            run("conda", "build", recipe, "--override-channels",
                "-c", "conda-forge", "-c", "bioconda", "--no-anaconda-upload",
                "--output-folder", channel, log=output / "conda-build.log",
                env={**os.environ, "CONDA_CHANNEL_PRIORITY": "strict"})
            packages = list(channel.glob("linux-64/metapathways-*.conda"))
            packages += list(channel.glob("linux-64/metapathways-*.tar.bz2"))
            if len(packages) != 1:
                raise ValueError(f"Expected one Linux package, found {packages}")
            shutil.copy2(packages[0], output / packages[0].name)
            run("conda", "index", channel, log=output / "conda-index.log")
            environment = work / "validation-env"
            run("conda", "create", "--yes", "--prefix", environment,
                "--override-channels", "--strict-channel-priority",
                "-c", channel.as_uri(), "-c", "conda-forge", "-c", "bioconda",
                f"metapathways={value}", log=output / "environment-create.log",
                env={**os.environ, "CONDA_ADD_PIP_AS_PYTHON_DEPENDENCY": "false"})
            runner = ["conda", "run", "--no-capture-output", "--prefix", environment]
            testdir = work / "integration"
            testdir.mkdir()
            run(*runner, "python", "-c",
                "from metapathways._version import __version__; "
                f"assert __version__ == {value!r}, __version__",
                cwd=testdir, log=output / "version-check.log")
            security_log = output / "security-dependencies.log"
            run(*runner, "python", "-c",
                "import sys, json, importlib.util; from importlib.metadata import version; "
                "print(json.dumps({'python': list(sys.version_info[:3]), "
                "'pip_present': importlib.util.find_spec('pip') is not None, "
                "'versions': {n: version(n) for n in ['urllib3', 'setuptools']}}))",
                cwd=testdir, log=security_log)
            validate_runtime_dependencies(json.loads(security_log.read_text()))
            run(*runner, "magsplitter", "--help", cwd=testdir, log=output / "magsplitter.log")
            run(*runner, "python", "-c", "import camelot_frs", cwd=testdir, log=output / "camelot.log")
            run(*runner, "metapathways", "prepare_test", "-o", "test", cwd=testdir, log=output / "test-inputs.log")
            run(*runner, "metapathways", "version", cwd=testdir,
                log=output / "cli-version.log")
            run(*runner, "python", snapshot / "scripts/check_installed_assets.py",
                cwd=testdir, log=output / "installed-assets.log")
            run(*runner, "metapathways", "build_db", "--test", "--memory", "2 GB", "--max_memory", "4 GB", "--max_cpus", "2", cwd=testdir,
                log=output / "build-db.log")
            run(*runner, "metapathways", "run", "--test", "--memory", "2 GB", "--max_memory", "4 GB", "--max_cpus", "2", cwd=testdir,
                log=output / "pipeline.log")
            validation.update(validate_run(testdir), scope="core-integration")
            for name in ["metapathways_steps_log.txt", "errors_warnings_log.txt"]:
                shutil.copy2(testdir / "test/k12_test" / name, output / name)
            shutil.copy2(testdir / "test/global_errors_warnings.txt", output)
            (output / "conda-explicit.txt").write_text(
                run("conda", "list", "--prefix", environment, "--explicit", "--sha256", capture=True) + "\n")
            (output / "pip-freeze.txt").write_text(
                run(*runner, "python", "-c",
                    "from importlib.metadata import distributions; "
                    "print('\\n'.join(sorted(f'{d.metadata[\"Name\"]}=={d.version}' for d in distributions())))",
                    capture=True) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    hashes = {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()}
    manifest = {"version": value, "build_number": number, "release_tag": tag, "commit": commit, "files": hashes}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    hashes["manifest.json"] = digest(output / "manifest.json")
    (output / "SHA256SUMS").write_text("".join(f"{sha}  {name}\n" for name, sha in sorted(hashes.items())))
    print(f"Built {value}: {output}")


def verify_artifacts(args):
    output = Path(args.output).resolve()
    manifest = json.loads((output / "manifest.json").read_text())
    ref = getattr(args, "ref", None)
    if ref:
        expected_version, expected_build = tag_version(ref)
        revision = f"refs/tags/{ref}"
        source = run("git", "show", f"{revision}:metapathways/_version.py", capture=True)
        if expected_build:
            recipe = run("git", "show", f"{revision}:conda_recipe/meta_template.yaml", capture=True)
            if recipe_build(recipe) != expected_build or manifest.get("build_number") != expected_build:
                raise ValueError("Build revision differs from release tag.")
        if version_from_text(source) != expected_version:
            raise ValueError("Release tag differs from its source version.")
        expected_commit = run("git", "rev-parse", f"{revision}^{{commit}}", capture=True)
    else:
        expected_version = version()
        expected_commit = run("git", "rev-parse", "HEAD", capture=True)
    if manifest["version"] != expected_version:
        raise ValueError("Artifact version differs from the requested source.")
    if manifest["commit"] != expected_commit:
        raise ValueError("Artifacts were built from a different commit.")
    if ref and manifest.get("release_tag", ref) != ref:
        raise ValueError("Artifact release tag differs from requested tag.")
    expected = set(manifest["files"]) | {"manifest.json", "SHA256SUMS"}
    if {p.name for p in output.iterdir()} != expected:
        raise ValueError("Artifact directory contains missing or unverified files.")
    if "validation.json" not in manifest["files"]:
        raise ValueError("Validation receipt is missing from manifest.")
    if any(p.is_symlink() or not p.is_file() for p in output.iterdir()):
        raise ValueError("Artifacts must be regular files, not links or directories.")
    for name, sha in manifest["files"].items():
        if Path(name).name != name or digest(output / name) != sha:
            raise ValueError(f"Artifact checksum mismatch: {name}")
    sums = "".join(f"{sha}  {name}\n" for name, sha in sorted(
        dict(manifest["files"], **{"manifest.json": digest(output / "manifest.json")}).items()))
    if (output / "SHA256SUMS").read_text() != sums:
        raise ValueError("SHA256SUMS differs from manifest.")
    validation = json.loads((output / "validation.json").read_text())
    if validation.get("scope") != "core-integration" or set(validation.get("successful_stages", [])) != STAGES:
        raise ValueError("Publishing requires a successful full Conda integration build.")
    if validation.get("commit") != manifest["commit"] or validation.get("version") != manifest["version"]:
        raise ValueError("Validation provenance differs from manifest.")
    print("Artifact checksums, revision, and integration receipt verified.")
    return output, manifest


def github_release_record(repo, tag):
    # Listing includes drafts, so interrupted first uploads can be resumed.
    raw = run("gh", "api", "--paginate", "--slurp",
              f"repos/{repo}/releases?per_page=100", capture=True)
    for page in json.loads(raw):
        for item in page:
            if item["tag_name"] == tag:
                return item
    return None


def github_asset_digest(repo, asset):
    recorded = asset.get("digest") or ""
    if recorded.startswith("sha256:"):
        return recorded.removeprefix("sha256:")
    # Older assets may not have an API digest. Download privately to hash them.
    result = subprocess.run(
        ["gh", "api", "-H", "Accept: application/octet-stream",
         f"repos/{repo}/releases/assets/{asset['id']}"],
        check=True, stdout=subprocess.PIPE)
    return hashlib.sha256(result.stdout).hexdigest()


def github_assets(output, value, directory):
    """Keep zero-byte logs inside a reproducible, complete release ZIP."""
    files = sorted(output.iterdir())
    bundle = Path(directory) / f"metapathways-{value}-release.zip"
    with zipfile.ZipFile(bundle, "w", compression=zipfile.ZIP_STORED) as archive:
        for path in files:
            info = zipfile.ZipInfo(path.name, date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes())
    # GitHub rejects zero-byte standalone release assets with Bad Content-Length.
    return [path for path in files if path.stat().st_size > 0] + [bundle]


def github_release(args):
    args.ref = args.tag
    output, manifest = verify_artifacts(args)
    repo = os.environ.get("GH_REPO", REPOSITORY)
    tag = args.tag
    record = github_release_record(repo, tag)
    if record is None:
        # Retain the create response: a newly created draft need not appear in
        # the collection listing immediately. All subsequent operations use ID.
        record = json.loads(run("gh", "api", "--method", "POST", f"repos/{repo}/releases",
            "-f", f"tag_name={tag}", "-f", f"name=MetaPathways {tag}",
            "-F", "draft=true", "-F", "generate_release_notes=true",
            "-F", f"prerelease={'true' if 'rc' in manifest['version'] else 'false'}",
            "-f", "body=The complete release ZIP includes all artifacts, validation logs "
            "(including empty success logs), and SHA256SUMS.", capture=True))
    release_id = record["id"]
    endpoint = f"repos/{repo}/releases/{release_id}"
    if record["tag_name"] != tag:
        raise ValueError("GitHub returned a release for a different tag.")
    with tempfile.TemporaryDirectory(prefix="metapathways-upload-") as temp:
        assets = github_assets(output, manifest["version"], temp)
        existing = {item["name"]: item for item in record["assets"]}
        for path in assets:
            if path.name in existing:
                if github_asset_digest(repo, existing[path.name]) != digest(path):
                    raise ValueError(f"Existing GitHub asset differs: {path.name}; refusing to overwrite.")
                print(f"Already uploaded and verified: {path.name}")
            else:
                upload = f"https://uploads.github.com/repos/{repo}/releases/{release_id}/assets"
                run("gh", "api", "--method", "POST", f"{upload}?{urlencode({'name': path.name})}",
                    "-H", "Content-Type: application/octet-stream", "--input", path, capture=True)
        record = json.loads(run("gh", "api", endpoint, capture=True))
        uploaded = {item["name"]: item for item in record["assets"]}
        for path in assets:
            if path.name not in uploaded or github_asset_digest(repo, uploaded[path.name]) != digest(path):
                raise ValueError(f"Uploaded GitHub asset failed checksum verification: {path.name}")
        if record["draft"]:
            run("gh", "api", "--method", "PATCH", endpoint, "-F", "draft=false", capture=True)
    print(f"Published and verified {tag}: https://github.com/{repo}/releases/tag/{tag}")


def publish(args):
    clean()
    value = version()
    branch = run("git", "branch", "--show-current", capture=True)
    target_branch = getattr(args, "branch", "main")
    if target_branch not in ("main", "dev") or branch != target_branch:
        raise ValueError(f"Publish from the {target_branch} branch after PR review and testing.")
    url = run("git", "remote", "get-url", "--push", args.remote, capture=True)
    if url.removesuffix(".git").rstrip("/") not in (
        f"git@github.com:{REPOSITORY}", f"https://github.com/{REPOSITORY}",
        f"ssh://git@github.com/{REPOSITORY}",
    ):
        raise ValueError(f"Remote must target {REPOSITORY}; got {url}")
    tag = release_tag(value, recipe_build((ROOT / "conda_recipe/meta_template.yaml").read_text()))
    if run("git", "ls-remote", args.remote, f"refs/tags/{tag}", capture=True):
        raise ValueError(f"{tag} already exists remotely. Rerun its CI job or prepare a new version.")
    existing = subprocess.run(["git", "show-ref", "--verify", "--quiet", f"refs/tags/{tag}"], cwd=ROOT)
    if existing.returncode == 0:
        if run("git", "rev-parse", f"{tag}^{{commit}}", capture=True) != run("git", "rev-parse", "HEAD", capture=True):
            raise ValueError(f"Local {tag} points to another commit; it will not be moved.")
        if run("git", "cat-file", "-t", tag, capture=True) != "tag":
            raise ValueError(f"Local {tag} must be an annotated tag.")
    else:
        run("git", "tag", "-a", tag, "-m", f"MetaPathways {value}")
    run("git", "push", "--atomic", args.remote, f"HEAD:refs/heads/{target_branch}", f"refs/tags/{tag}")
    print(f"CI will build and test {tag}; publication requires manual selection: https://github.com/{REPOSITORY}/actions")


def upload_conda(args):
    output, manifest = verify_artifacts(args)
    packages = [output / n for n in manifest["files"] if n.endswith((".conda", ".tar.bz2"))]
    if len(packages) != 1:
        raise ValueError("Expected exactly one validated Conda package.")
    label = "rc" if "rc" in manifest["version"] else "main"
    # anaconda-client reads BINSTAR_API_TOKEN; never put credentials on the command line.
    run("anaconda", "upload", "--user", os.environ.get("ANACONDA_OWNER", "hallamlab"), "--label", label, packages[0])


def container_command(args):
    import containers
    getattr(containers, args.container_action)(args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare", help="Set an explicit stable or release-candidate version.")
    p.add_argument("version")
    p.add_argument("--build-number", type=int,
                   help="Default: preserve the current build for the same version; 0 for a new version.")
    p.set_defaults(func=prepare)
    p = commands.add_parser("build", help="Build archives and a Conda package; test the installed package.")
    p.add_argument("--output", default=str(ROOT / "dist/release"))
    p.add_argument("--source-only", action="store_true", help="Build source archives only; not publishable.")
    p.add_argument("--tag", help="CI: require this tag to match the version and checked-out commit.")
    p.set_defaults(func=build)
    for name, func in [("verify-artifacts", verify_artifacts), ("upload-conda", upload_conda)]:
        p = commands.add_parser(name)
        p.add_argument("--output", default=str(ROOT / "dist/release"))
        p.add_argument("--ref", help="Verify against this existing release tag instead of HEAD.")
        p.set_defaults(func=func)
    p = commands.add_parser("github-release", help="Publish or resume a release; preserve empty logs in a ZIP.")
    p.add_argument("tag")
    p.add_argument("--output", default=str(ROOT / "dist/release"))
    p.set_defaults(func=github_release)
    p = commands.add_parser("publish", help="Push the reviewed main/dev branch and its release tag using existing Git credentials.")
    p.add_argument("--remote", default="origin")
    p.add_argument("--branch", choices=("main", "dev"), default="main",
                   help="Reviewed release branch [main]; never a feature branch.")
    p.set_defaults(func=publish)
    for command, action in [("container-build", "build"), ("container-push", "push"),
                            ("container-verify", "verify"), ("container-attach", "attach"), ("quay-description", "description")]:
        p = commands.add_parser(command)
        p.add_argument("--output", default=str(ROOT / "dist/release"))
        p.add_argument("--container-output", default=str(ROOT / "dist/containers"))
        p.add_argument("--ref", help="Verify against an existing release tag.")
        p.set_defaults(func=container_command, container_action=action)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Release stopped: {error}\n")


if __name__ == "__main__":
    main()
