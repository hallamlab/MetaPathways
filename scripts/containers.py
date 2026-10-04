#!/usr/bin/env python3
"""Build tested Docker/Apptainer releases and update the Quay repository overview."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request
from urllib.parse import urlsplit

import release

IMAGE = "quay.io/hallamlab/metapathways"
QUAY_API = "https://quay.io/api/v1/repository/hallamlab/metapathways"


def explicit_lock(text, package):
    """Preserve the tested solve; replace only the local release package URL."""
    lines = text.splitlines()
    if "@EXPLICIT" not in lines:
        raise ValueError("Missing explicit Conda dependency export")
    matches = 0
    for i, line in enumerate(lines):
        if not line or line.startswith(("#", "@")):
            continue
        url = urlsplit(line)
        if Path(url.path).name == package:
            matches += 1
            lines[i] = f"file:///tmp/release-package/{package}"
        elif url.scheme != "https" or url.hostname != "conda.anaconda.org":
            raise ValueError("Unexpected dependency URL in explicit export")
    if matches != 1:
        raise ValueError("Explicit export must contain exactly one release package")
    return "\n".join(lines) + "\n"


def build(args):
    source, manifest = release.verify_artifacts(args)
    output = Path(args.container_output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError(f"Container output must be empty: {output}")
    packages = [n for n in manifest["files"] if n.endswith((".conda", ".tar.bz2"))]
    if len(packages) != 1:
        raise ValueError("Expected one validated Conda package")
    package = packages[0]
    value = manifest["version"]
    image = f"metapathways-release:{manifest['commit'][:12]}"
    with tempfile.TemporaryDirectory(prefix="metapathways-container-") as temp:
        work = Path(temp)
        context = work / "context"
        (context / "package").mkdir(parents=True)
        shutil.copy2(release.ROOT / "docker/Dockerfile.release", context / "Dockerfile")
        shutil.copy2(source / package, context / "package" / package)
        (context / "package/explicit.txt").write_text(
            explicit_lock((source / "conda-explicit.txt").read_text(), package))
        release.run("docker", "build", "--pull", "--no-cache", "--platform", "linux/amd64", "--build-arg", f"VERSION={value}",
                    "--build-arg", f"REVISION={manifest['commit']}", "--tag", image, context,
                    log=output / "docker-build.log")
        integration = work / "integration"
        integration.mkdir(mode=0o777)
        integration.chmod(0o777)
        release.run("docker", "run", "--rm", "--volume", f"{integration}:/work", image,
                    "bash", "-c", "umask 0000; trap 'find /work -mindepth 1 -exec chmod a+rwX {} +' EXIT; "
                    "metapathways build_db --test --memory '2 GB' --max_memory '4 GB' --max_cpus 2 && metapathways run --test --memory '2 GB' --max_memory '4 GB' --max_cpus 2",
                    log=output / "docker-integration.log")
        receipt = release.validate_run(integration)
        for name in ["metapathways_steps_log.txt", "errors_warnings_log.txt"]:
            shutil.copy2(integration / "test/k12_test" / name, output / name)
        # Docker archive conversion uses the very same locally tested image.
        archive = work / "image.tar"
        release.run("docker", "save", "--output", archive, image)
        sif = output / f"metapathways-{value}-linux-amd64.sif"
        release.run("apptainer", "build", sif, f"docker-archive://{archive}",
                    log=output / "apptainer-build.log")
        release.run("apptainer", "exec", "--cleanenv", sif, "python", "-c",
                    "from metapathways._version import __version__; "
                    f"assert __version__ == {value!r}, __version__",
                    log=output / "apptainer-version.log")
        release.run("apptainer", "exec", "--cleanenv", sif, "metapathways", "run", "--help",
                    log=output / "apptainer-cli.log")
        receipt.update(version=value, release_tag=manifest.get("release_tag", f"v{value}"),
                       commit=manifest["commit"], package=package,
                       package_sha256=release.digest(source / package), image=image,
                       image_id=release.run("docker", "image", "inspect", "--format", "{{.Id}}", image, capture=True),
                       sif=sif.name, sif_sha256=release.digest(sif),
                       apptainer_validation="version assertion and run --help; Docker full integration")
        (output / "container-validation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    (output / "container-SHA256SUMS").write_text("".join(
        f"{release.digest(p)}  {p.name}\n" for p in sorted(output.iterdir()) if p.is_file()))
    print(f"Validated Docker image and Apptainer SIF: {output}")


def verify(args):
    _, manifest = release.verify_artifacts(args)
    output = Path(args.container_output).resolve()
    sums = {}
    for line in (output / "container-SHA256SUMS").read_text().splitlines():
        sha, name = line.split("  ", 1)
        if Path(name).name != name or name in sums:
            raise ValueError("Invalid container checksum entry")
        sums[name] = sha
    if set(p.name for p in output.iterdir()) != set(sums) | {"container-SHA256SUMS"}:
        raise ValueError("Unexpected container artifact files")
    for name, sha in sums.items():
        path = output / name
        if path.is_symlink() or release.digest(path) != sha:
            raise ValueError(f"Container checksum mismatch: {name}")
    receipt = json.loads((output / "container-validation.json").read_text())
    if (receipt["commit"], receipt["version"]) != (manifest["commit"], manifest["version"]):
        raise ValueError("Container source differs from release")
    if receipt.get("release_tag", f"v{receipt['version']}") != manifest.get("release_tag", f"v{manifest['version']}"):
        raise ValueError("Container release tag differs from package release")
    if manifest["files"].get(receipt["package"]) != receipt["package_sha256"]:
        raise ValueError("Container package differs from release")
    if set(receipt["successful_stages"]) != release.STAGES or set(receipt["nonempty_outputs"]) != set(release.OUTPUTS):
        raise ValueError("Container integration validation is incomplete")
    if sums.get(receipt["sif"]) != receipt["sif_sha256"]:
        raise ValueError("SIF validation checksum mismatch")
    return output, receipt


def push(args):
    _, receipt = verify(args)
    image = receipt["image"]
    if release.run("docker", "image", "inspect", "--format", "{{.Id}}", image, capture=True) != receipt["image_id"]:
        raise ValueError("Local image changed since validation")
    value = receipt["version"]
    tag = receipt.get("release_tag", f"v{value}")
    tags = list(dict.fromkeys([tag.removeprefix("v"), tag, value, f"v{value}"]))
    tags += [] if "rc" in value else ["latest"]
    for tag in tags:
        target = f"{os.environ.get('QUAY_REPOSITORY', IMAGE)}:{tag}"
        release.run("docker", "tag", image, target)
        release.run("docker", "push", target)


def attach(args):
    output, receipt = verify(args)
    repo = os.environ.get("GH_REPO", release.REPOSITORY)
    tag = receipt.get("release_tag", f"v{receipt['version']}")
    record = release.github_release_record(repo, tag)
    if record is None or record["draft"]:
        raise ValueError("Publish the GitHub release first")
    existing = {a["name"]: a for a in record["assets"]}
    # Bundle logs so empty success logs are never standalone GitHub uploads.
    with tempfile.TemporaryDirectory(prefix="metapathways-container-upload-") as temp:
        bundle = Path(temp) / f"metapathways-{receipt['version']}-container-validation.zip"
        with release.zipfile.ZipFile(bundle, "w", release.zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(output.iterdir()):
                if path.suffix != ".sif" and path.name != "docker-image.tar.gz":
                    info = release.zipfile.ZipInfo(path.name, (1980, 1, 1, 0, 0, 0))
                    archive.writestr(info, path.read_bytes())
        for path in [output / receipt["sif"], output / "container-SHA256SUMS", bundle]:
            if path.name in existing:
                if release.github_asset_digest(repo, existing[path.name]) != release.digest(path):
                    raise ValueError(f"Existing container asset differs: {path.name}; refusing overwrite")
            else:
                release.run("gh", "release", "upload", tag, path, "--repo", repo)
        uploaded = {a["name"]: a for a in release.github_release_record(repo, tag)["assets"]}
        for path in [output / receipt["sif"], output / "container-SHA256SUMS", bundle]:
            if path.name not in uploaded or release.github_asset_digest(repo, uploaded[path.name]) != release.digest(path):
                raise ValueError(f"Container upload verification failed: {path.name}")


def description(args):
    token = os.environ.get("QUAY_API_TOKEN", "")
    if not token:
        raise ValueError("Set QUAY_API_TOKEN to a Quay OAuth token with repo:write access")
    text = (release.ROOT / "docker/README.quay.md").read_text()
    request = urllib.request.Request(QUAY_API, method="PUT",
        data=json.dumps({"description": text}).encode(),
        headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        response.read()
    with urllib.request.urlopen(QUAY_API, timeout=60) as response:
        if json.load(response)["description"] != text:
            raise ValueError("Quay overview verification failed")
    print("Updated and verified https://quay.io/repository/hallamlab/metapathways")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, func in [("build", build), ("push", push), ("attach", attach), ("verify", verify)]:
        p = commands.add_parser(name)
        p.add_argument("--output", default=str(release.ROOT / "dist/release"))
        p.add_argument("--container-output", default=str(release.ROOT / "dist/containers"))
        p.add_argument("--ref", help="Existing release tag for artifact verification")
        p.set_defaults(func=func)
    commands.add_parser("description").set_defaults(func=description)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Container release stopped: {error}\n")


if __name__ == "__main__":
    main()
