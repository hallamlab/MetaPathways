#!/bin/bash
#SBATCH --job-name=mp-controller
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=mp-controller-%j.log

set -euo pipefail
# Submit with: sbatch --account=YOUR_ACCOUNT analysis-controller.sbatch run-analysis.sh
# Activate the MetaPathways environment in run-analysis.sh and put the complete
# metapathways analysis_wf command there, including --executor slurm --account.
# sbatch, squeue and scancel must be usable from this compute node.
# This limits the Java heap, not total JVM memory or native thread counts.
export NXF_OPTS="${NXF_OPTS:--Xmx4g}"
if [[ $# -ne 1 ]]; then
    echo 'Usage: sbatch --account=ACCOUNT analysis-controller.sbatch run-analysis.sh' >&2
    exit 2
fi
bash "$1"
