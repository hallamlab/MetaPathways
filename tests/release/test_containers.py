"""Container provenance and dependency-lock safety checks."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import containers


class ContainerTests(unittest.TestCase):
    def test_lock_preserves_solve_and_replaces_only_release_package(self):
        text = "# explicit\n@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/python.conda#abc\nfile:///old/channel/metapathways.tar.bz2\n"
        result = containers.explicit_lock(text, "metapathways.tar.bz2")
        self.assertIn("python.conda#abc", result)
        self.assertIn("file:///tmp/release-package/metapathways.tar.bz2", result)
        for invalid in [text.replace("@EXPLICIT", ""), text + "file:///other/package.conda\n",
                        text + "file:///duplicate/metapathways.tar.bz2\n",
                        text.replace("conda.anaconda.org", "untrusted.invalid")]:
            with self.assertRaises(ValueError):
                containers.explicit_lock(invalid, "metapathways.tar.bz2")

    def test_verify_rejects_tampered_sif_and_wrong_source(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "image.sif").write_bytes(b"test SIF")
            receipt = dict(commit="abc", version="3.5.0", package="package.conda",
                           package_sha256="package-sha", sif="image.sif",
                           sif_sha256=containers.release.digest(root / "image.sif"),
                           successful_stages=sorted(containers.release.STAGES),
                           nonempty_outputs=containers.release.OUTPUTS)
            (root / "container-validation.json").write_text(json.dumps(receipt))
            (root / "container-SHA256SUMS").write_text("".join(
                f"{containers.release.digest(p)}  {p.name}\n" for p in sorted(root.iterdir())))
            manifest = dict(commit="abc", version="3.5.0", files={"package.conda": "package-sha"})
            args = SimpleNamespace(container_output=temp)
            with patch.object(containers.release, "verify_artifacts", return_value=(root, manifest)):
                containers.verify(args)
                manifest["commit"] = "wrong"
                with self.assertRaisesRegex(ValueError, "source differs"):
                    containers.verify(args)
                manifest["commit"] = "abc"
                (root / "image.sif").write_bytes(b"tampered")
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    containers.verify(args)

    def test_description_requires_separate_api_credential(self):
        with patch.dict(containers.os.environ, {}, clear=True):
            with self.assertRaisesRegex(ValueError, "QUAY_API_TOKEN"):
                containers.description(SimpleNamespace())
