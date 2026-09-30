"""Release safety checks; no credentials or remote services are used."""
import importlib.util
import json
import subprocess
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/release.py"
spec = importlib.util.spec_from_file_location("release", SCRIPT)
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class ReleaseTests(unittest.TestCase):
    def test_explicit_versions_only(self):
        for value in ["3.5.0", "v3.5.1", "3.6.0rc1"]:
            self.assertEqual(release.valid_version(value), value.removeprefix("v"))
        for value in ["3.5.0.dev64", "3.5", "../3.5.0", "3.5.0;echo bad", "03.5.0", "3.5.0rc0"]:
            with self.assertRaises(ValueError):
                release.valid_version(value)

    def test_preparation_preserves_unrelated_changes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "metapathways").mkdir()
            (root / "conda_recipe").mkdir()
            (root / "metapathways/_version.py").write_text(
                '__version__ = "3.5.0"\n__status__ = "Release"\n')
            (root / "conda_recipe/meta_template.yaml").write_text("build:\n  number: 0\n")
            (root / "README.md").write_text(
                "User edit\nhttps://img.shields.io/badge/Version-3.5-blue.svg)\n")
            with patch.object(release, "ROOT", root):
                release.prepare(SimpleNamespace(version="3.5.1rc1", build_number=2))
            self.assertEqual(release.version(root), "3.5.1rc1")
            self.assertIn("User edit", (root / "README.md").read_text())
            self.assertIn("number: 2", (root / "conda_recipe/meta_template.yaml").read_text())
            with patch.object(release, "ROOT", root):
                release.prepare(SimpleNamespace(version="3.5.1rc1", build_number=None))
            self.assertIn("number: 2", (root / "conda_recipe/meta_template.yaml").read_text())
            with patch.object(release, "ROOT", root):
                release.prepare(SimpleNamespace(version="3.5.1", build_number=None))
            self.assertIn("number: 0", (root / "conda_recipe/meta_template.yaml").read_text())

    def fixture(self, root):
        sample = root / "test/k12_test"
        sample.mkdir(parents=True)
        (sample / "metapathways_steps_log.txt").write_text(
            "".join(f"{stage}\tSUCCESS - Time elapsed: 1\n" for stage in release.STAGES))
        (sample / "errors_warnings_log.txt").write_text("#STEP\tRPKM_CALCULATION\n")
        (root / "test/global_errors_warnings.txt").write_text("")
        for name in release.OUTPUTS:
            path = sample / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("data\n")
        return sample

    def test_integration_checks_content_not_exit_status(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sample = self.fixture(root)
            self.assertEqual(len(release.validate_run(root)["successful_stages"]), 16)
            output = sample / release.OUTPUTS[0]
            output.write_text("")
            with self.assertRaisesRegex(ValueError, "empty output"):
                release.validate_run(root)
            output.write_text("data")
            (sample / "errors_warnings_log.txt").write_text("ERROR\ttool failed\n")
            with self.assertRaisesRegex(ValueError, "Review errors"):
                release.validate_run(root)

    def test_duplicate_success_cannot_hide_missing_stage(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sample = self.fixture(root)
            path = sample / "metapathways_steps_log.txt"
            path.write_text(path.read_text().replace("COMPUTE_TPM\t", "ORF_PREDICTION\t"))
            with self.assertRaisesRegex(ValueError, "incomplete"):
                release.validate_run(root)

    def test_dirty_checkout_stops_build_before_artifacts(self):
        with patch.object(release, "run", return_value=" M README.md"):
            with self.assertRaisesRegex(ValueError, "Commit or stash"):
                release.build(SimpleNamespace())

    def test_existing_remote_tag_never_overwritten(self):
        responses = ["", "dev", "git@github.com:hallamlab/MetaPathways.git",
                     "abc refs/tags/v3.5.0"]
        with patch.object(release, "version", return_value="3.5.0"), patch.object(
            release, "run", side_effect=responses
        ) as run:
            with self.assertRaisesRegex(ValueError, "already exists remotely"):
                release.publish(SimpleNamespace(remote="origin"))
            self.assertFalse(any("push" in c.args or "tag" in c.args for c in run.call_args_list))

    def test_publish_pushes_annotated_tag_and_branch_atomically(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            remote, checkout = root / "remote.git", root / "checkout"
            subprocess.run(["git", "init", "--bare", str(remote)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run(["git", "init", "-b", "dev", str(checkout)], check=True,
                           stdout=subprocess.DEVNULL)
            (checkout / "README.md").write_text("Release fixture\n")
            for command in [
                ["git", "add", "."],
                ["git", "-c", "user.name=Release Test", "-c", "user.email=test@example.invalid",
                 "commit", "-qm", "Fixture"],
                ["git", "config", "user.name", "Release Test"],
                ["git", "config", "user.email", "test@example.invalid"],
                ["git", "remote", "add", "origin", str(remote)],
            ]:
                subprocess.run(command, cwd=checkout, check=True)
            original_run = release.run

            def local_run(*args, **kwargs):
                if args[:3] == ("git", "remote", "get-url"):
                    return "git@github.com:hallamlab/MetaPathways.git"
                return original_run(*args, cwd=checkout, **kwargs)

            with patch.object(release, "ROOT", checkout), patch.object(
                release, "version", return_value="3.5.0"
            ), patch.object(release, "run", side_effect=local_run):
                release.publish(SimpleNamespace(remote="origin"))
                with self.assertRaisesRegex(ValueError, "already exists remotely"):
                    release.publish(SimpleNamespace(remote="origin"))
            kind = subprocess.check_output(
                ["git", "cat-file", "-t", "v3.5.0"], cwd=remote, text=True).strip()
            self.assertEqual(kind, "tag")
            branch = subprocess.check_output(
                ["git", "rev-parse", "dev"], cwd=remote, text=True).strip()
            tagged = subprocess.check_output(
                ["git", "rev-parse", "v3.5.0^{commit}"], cwd=remote, text=True).strip()
            self.assertEqual(branch, tagged)

    def test_manifest_rejects_tampering_and_extra_files(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            receipt = {"scope": "core-integration", "version": "3.5.0", "commit": "abc",
                       "successful_stages": sorted(release.STAGES)}
            path = root / "validation.json"
            path.write_text(json.dumps(receipt))
            manifest = {"version": "3.5.0", "commit": "abc",
                        "files": {"validation.json": release.digest(path)}}
            (root / "manifest.json").write_text(json.dumps(manifest))
            sums = dict(manifest["files"], **{"manifest.json": release.digest(root / "manifest.json")})
            (root / "SHA256SUMS").write_text(
                "".join(f"{sha}  {name}\n" for name, sha in sorted(sums.items())))
            with patch.object(release, "version", return_value="3.5.0"), patch.object(
                release, "run", return_value="abc"
            ):
                release.verify_artifacts(SimpleNamespace(output=root))
                extra = root / "unverified.txt"
                extra.write_text("Do not upload this")
                with self.assertRaisesRegex(ValueError, "unverified files"):
                    release.verify_artifacts(SimpleNamespace(output=root))
                extra.unlink()
                path.write_text("{}")
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    release.verify_artifacts(SimpleNamespace(output=root))
                receipt["scope"] = "source-only"
                path.write_text(json.dumps(receipt))
                manifest["files"]["validation.json"] = release.digest(path)
                (root / "manifest.json").write_text(json.dumps(manifest))
                sums = dict(manifest["files"], **{"manifest.json": release.digest(root / "manifest.json")})
                (root / "SHA256SUMS").write_text(
                    "".join(f"{sha}  {name}\n" for name, sha in sorted(sums.items())))
                with self.assertRaisesRegex(ValueError, "successful full Conda"):
                    release.verify_artifacts(SimpleNamespace(output=root))


if __name__ == "__main__":
    unittest.main()
