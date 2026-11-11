# Maintainer CLI reference

[Development](development.md) · [Reaction screening](ptools-maintenance.md)

Run `metapathways COMMAND --help` for the help of your installed revision.

## screen_pt

```text
usage: metapathways [-h] -d REFDB_DIR -o OUTPUT_DIR [--image IMAGE] [--publish]
                    [--reactions REACTIONS [REACTIONS ...]] [--batch_size BATCH_SIZE]
                    [--max_tasks MAX_TASKS] [--confirm_runs CONFIRM_RUNS] [--timeout TIMEOUT]
                    [--scratch_dir SCRATCH_DIR]

Screen explicit MPDB reaction assignments against a licensed PTools SIF

options:
  -h, --help            show this help message and exit
  -d REFDB_DIR, --refdb_dir REFDB_DIR
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
  --image IMAGE         PTools SIF [registered build_pt image]
  --publish             Save a completed full screen compatibility list in the MPDB
  --reactions REACTIONS [REACTIONS ...]
                        Optional reaction IDs for a targeted screen
  --batch_size BATCH_SIZE
  --max_tasks MAX_TASKS
                        Concurrent isolated containers [1]
  --confirm_runs CONFIRM_RUNS
  --timeout TIMEOUT     Seconds per build [1800]; timeout is inconclusive
  --scratch_dir SCRATCH_DIR
                        Local temporary storage for private PGDB builds
```
