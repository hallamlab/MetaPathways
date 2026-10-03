## Problem and resulting behavior

Describe the concrete change and any compatibility implications.

## Validation

- Tested commit:
- Local tests and installed-package checks:
- Single/two-sample reviewer workflow and checkpoint reuse:
- Report/explorer and CSV export:
- Optional Pathway Tools and Slurm checks (or explicitly untested):
- Nonpublishing Release workflow/artifact links:

See [the tester checklist](https://metapathways.readthedocs.io/en/latest/pr-testing.html). Feature PRs target `dev`; keep them in draft until single-server and Slurm benchmark review plus independent tester sign-off are complete. Production promotion uses a separate `dev` → `master` PR. Do not include licensed installers/images, credentials, or production outputs.
