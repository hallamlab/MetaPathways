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

    def test_build_revision_keeps_application_version(self):
        self.assertEqual(release.release_tag("3.5.1", 1), "v3.5.1-build1")
        self.assertEqual(release.tag_version("v3.5.1-build1"), ("3.5.1", 1))
        self.assertEqual(release.tag_version("v3.5.1"), ("3.5.1", 0))
        for bad in ["3.5.1", "v3.5.1-build0", "v3.5.1-build01", "v3.5.1-build-1"]:
            with self.assertRaises(ValueError):
                release.tag_version(bad)

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
            (root / "CITATION.cff").write_text('title: MetaPathways\nversion: "3.5.0"\n')
            with patch.object(release, "ROOT", root):
                release.prepare(SimpleNamespace(version="3.5.1rc1", build_number=2))
            self.assertEqual(release.version(root), "3.5.1rc1")
            self.assertIn('version: "3.5.1rc1"', (root / "CITATION.cff").read_text())
            self.assertIn("User edit", (root / "README.md").read_text())
            self.assertIn("number: 2", (root / "conda_recipe/meta_template.yaml").read_text())
            with patch.object(release, "ROOT", root):
                release.prepare(SimpleNamespace(version="3.5.1rc1", build_number=None))
            self.assertIn("number: 2", (root / "conda_recipe/meta_template.yaml").read_text())
            with patch.object(release, "ROOT", root):
                release.prepare(SimpleNamespace(version="3.5.1", build_number=None))
            self.assertIn("number: 0", (root / "conda_recipe/meta_template.yaml").read_text())
            self.assertEqual((root / "CITATION.cff").read_text(), 'title: MetaPathways\nversion: "3.5.1"\n')

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
        with patch.object(release, "recipe_build", return_value=0), patch.object(release, "version", return_value="3.5.0"), patch.object(
            release, "run", side_effect=responses
        ) as run:
            with self.assertRaisesRegex(ValueError, "already exists remotely"):
                release.publish(SimpleNamespace(remote="origin"))
            self.assertFalse(any("push" in c.args or "tag" in c.args for c in run.call_args_list))

    def test_feature_branch_cannot_publish(self):
        with patch.object(release, "version", return_value="3.5.2"), patch.object(
            release, "run", side_effect=["", "feat/nextflow-controller-db-build"]
        ) as run:
            with self.assertRaisesRegex(ValueError, "dev branch"):
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
            (checkout / "conda_recipe").mkdir()
            (checkout / "conda_recipe/meta_template.yaml").write_text("build:\n  number: 0\n")
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
            with patch.object(release, "recipe_build", return_value=0), patch.object(release, "version", return_value="3.5.0"), patch.object(
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

    def test_github_bundle_preserves_empty_logs_and_is_repeatable(self):
        import zipfile
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            output = root / "artifacts"
            output.mkdir()
            (output / "global_errors_warnings.txt").write_bytes(b"")
            (output / "package.tar.gz").write_bytes(b"package")
            first = release.github_assets(output, "3.5.0", root)
            first_digest = release.digest(first[-1])
            self.assertTrue(all(p.stat().st_size > 0 for p in first))
            self.assertNotIn("global_errors_warnings.txt", [p.name for p in first])
            with zipfile.ZipFile(first[-1]) as archive:
                self.assertEqual(archive.read("global_errors_warnings.txt"), b"")
                self.assertEqual(archive.read("package.tar.gz"), b"package")
            (output / "package.tar.gz").touch()
            second = release.github_assets(output, "3.5.0", root)
            self.assertEqual(first_digest, release.digest(second[-1]))

    def test_github_recovery_resumes_partial_draft_without_overwrites(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "global_errors_warnings.txt").write_bytes(b"")
            package = root / "package.tar.gz"
            package.write_bytes(b"package")
            record = {"id": 123, "tag_name": "v3.5.0", "draft": True, "assets": [
                {"name": package.name, "digest": "sha256:" + release.digest(package)}
            ]}
            uploaded = []

            def fake_run(*args, **kwargs):
                if args[:4] == ("gh", "api", "--method", "POST"):
                    path = Path(args[args.index("--input") + 1])
                    self.assertGreater(path.stat().st_size, 0)
                    self.assertNotIn("--clobber", args)
                    uploaded.append(path.name)
                    record["assets"].append({"name": path.name,
                                             "digest": "sha256:" + release.digest(path)})
                elif args[:4] == ("gh", "api", "--method", "PATCH"):
                    record["draft"] = False
                return json.dumps(record)

            with patch.object(release, "verify_artifacts", return_value=(root, {"version": "3.5.0"})), patch.object(
                release, "github_release_record", return_value=record
            ), patch.object(release, "run", side_effect=fake_run):
                release.github_release(SimpleNamespace(tag="v3.5.0", output=root))
                self.assertEqual(uploaded, ["metapathways-3.5.0-release.zip"])
                self.assertFalse(record["draft"])
                # A completed release can be rerun without uploading anything again.
                release.github_release(SimpleNamespace(tag="v3.5.0", output=root))
                self.assertEqual(len(uploaded), 1)
                record["assets"][0]["digest"] = "sha256:wrong"
                with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                    release.github_release(SimpleNamespace(tag="v3.5.0", output=root))

    def test_new_draft_uses_creation_id_without_relisting(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "package.tar.gz").write_bytes(b"package")
            record = {"id": 456, "tag_name": "v3.5.1", "draft": True, "assets": []}
            def fake_run(*args, **kwargs):
                if "--input" in args:
                    self.assertIn("/releases/456/assets?", args[4])
                    path = Path(args[args.index("--input") + 1])
                    record["assets"].append({"name": path.name,
                        "digest": "sha256:" + release.digest(path)})
                elif "PATCH" in args:
                    self.assertEqual(args[4], "repos/hallamlab/MetaPathways/releases/456")
                    record["draft"] = False
                return json.dumps(record)
            with patch.object(release, "verify_artifacts", return_value=(root, {"version": "3.5.1"})), \
                 patch.object(release, "github_release_record", return_value=None) as lookup, \
                 patch.object(release, "run", side_effect=fake_run):
                release.github_release(SimpleNamespace(tag="v3.5.1", output=root))
                lookup.assert_called_once()
                self.assertFalse(record["draft"])
                self.assertEqual(len(record["assets"]), 2)

    def test_recovery_checks_original_tag_commit_not_workflow_commit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {"version": "3.5.0", "commit": "original", "files": {}}
            (root / "manifest.json").write_text(json.dumps(manifest))
            with patch.object(release, "run", side_effect=[
                '__version__ = "3.5.0"', "other-commit"
            ]):
                with self.assertRaisesRegex(ValueError, "different commit"):
                    release.verify_artifacts(SimpleNamespace(output=root, ref="v3.5.0"))


if __name__ == "__main__":
    unittest.main()
