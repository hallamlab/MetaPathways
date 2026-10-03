# PR tester checklist

[Home](index.md) · [Reviewer commands](reviewer-test.md) · [Release process](releasing.md)

Test the feature branch before merging to `dev` or `master`. Use a clean installation and a fresh output directory. Record `git rev-parse HEAD`, the platform, and the installed dependency versions. A version number alone does not identify the tested commit.

1. Follow [source installation](getting-started.md). Check out the PR's exact commit after cloning; do not test an older editable installation by accident. Record the MAGSplitter and Camelot revisions installed separately.
2. Follow the [CAMI II reviewer walkthrough](reviewer-test.md) ([Meyer et al., 2022](#cami-references)), running its single-sample and two-sample commands with `--skip_ptools`. Check required task outcomes, nonempty abundance outputs, the three genome-bin assignments per sample, and distinct paired read paths.
3. Repeat the two-sample command unchanged. All successful workflow tasks, including `COMPUTE_TPM`, should report `ALREADY_COMPUTED`. Report rebuilding may still run. A repeat is not a fresh performance measurement.
4. Run all three inputs through automatic discovery with a new output directory. Confirm the same sample IDs and input associations as `all.tsv`. Test an intentionally mismatched filename in a separate copy: planning must stop with an actionable input-layout error.
5. Serve the report and open `EDA_portal.html`. Filter each sample, search annotations, and export a CSV. Confirm sample/entity IDs and absence of mixed sample records. Sparse annotations and empty RNA/pathway tables may be expected with the tiny reference fixture; skipped pathways must not appear as successful inference.
6. If licensed, build a Pathway Tools image from your installer. Repeat a small workflow into a fresh output with `--taxprune --taxonomic_scope all` and without `--skip_ptools`. Confirm community/MAG status, sequence-backed inputs, archives and pathway exports. Keep licensed inputs/images private. Record expected no-pathway outcomes separately from failures.
7. If a cluster is available, test a small submission using the documented Slurm settings. Confirm scheduler resource requests and bounded concurrent submissions. Otherwise explicitly mark Slurm untested.
8. Have a maintainer manually run the Release workflow on this branch with **blank `release_tag` and `source_run_id`**, optionally enabling `test_containers`. Review the Conda installed-package integration, container validation, and security reports. No release tags or public uploads are needed for this step.

Attach the tested commit, commands, task summaries, dependency versions, and relevant logs to the PR, with local credentials and private input paths removed as needed. State which optional checks were not performed. A maintainer should resolve failures and obtain sign-off on the final commit before merging; subsequent code changes require retesting the affected behavior.

```{include} includes/cami-references.md
```
