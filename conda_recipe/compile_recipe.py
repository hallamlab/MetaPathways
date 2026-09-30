"""Render the Conda recipe from the exact source archive and runtime environment."""
import argparse
import hashlib
from pathlib import Path
import shlex
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from setup import NAME, VERSION, ENTRY_POINTS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdist", type=Path, default=ROOT / "dist" / f"{NAME}-{VERSION}.tar.gz")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "conda_recipe")
    args = parser.parse_args()
    archive = args.sdist.resolve()
    if archive.name != f"{NAME}-{VERSION}.tar.gz" or not archive.is_file():
        parser.error(f"Expected source archive {NAME}-{VERSION}.tar.gz; got {archive}")
    deps = yaml.safe_load((ROOT / "docker/conda_base.yml").read_text())
    if not all(isinstance(d, str) for d in deps["dependencies"]):
        parser.error("Conda runtime dependencies must be explicit package strings, not pip/VCS entries.")
    replacements = {
        "<NAME>": NAME,
        "<VERSION>": VERSION,
        "<ENTRY>": "\n".join(f"    - {e}" for e in ENTRY_POINTS),
        "<REQUIREMENTS>": "\n".join(f"    - {d}" for d in deps["dependencies"]),
        "<TAR>": archive.as_uri(),
        "<SHA256>": hashlib.sha256(archive.read_bytes()).hexdigest(),
    }
    text = (ROOT / "conda_recipe/meta_template.yaml").read_text()
    for key, value in replacements.items():
        text = text.replace(key, value)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "meta.yaml").write_text(text)
    command = ["conda", "build", "--override-channels", "-c", "conda-forge", "-c", "bioconda",
               "--no-anaconda-upload", "--output-folder", str(ROOT / "conda_build"),
               str(args.output_dir.resolve())]
    wrapper = args.output_dir / "call_build.sh"
    wrapper.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + shlex.join(command) + "\n")
    wrapper.chmod(0o755)


if __name__ == "__main__":
    main()
