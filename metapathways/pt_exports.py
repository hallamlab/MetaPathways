"""Handle Pathway Tools' empty-pathway export without altering the raw PGDB."""
from pathlib import Path


EMPTY_EXPORT_ERROR = (
    "Error: The flat-file generation program didn't specify what kind of data "
    "to put in this file."
)
PATHWAY_COLUMNS = (
    'SAMPLE', 'PWY_NAME', 'PWY_COMMON_NAME', 'PWY_SCORE', 'NUM_REACTIONS',
    'NUM_COVERED_REACTIONS', 'ORF_COUNT', 'ORFS',
)
REPORT_COLUMNS = (
    'Pathway Name', 'Pathway Frame-id', 'Pathway Class Name',
    'Pathway Class Frame-id', 'Pathway Score', 'Pathway Frequency Score',
    'Pathway Abundance', 'Reason to Keep', 'Pathway URL',
)


def write_verified_empty_pathways(flatpath, outfile):
    """Return True only for an independently verified empty pathway export.

    Some Pathway Tools versions export an error sentence for an empty frame
    class. Do not pass this sentence to Camelot or silently discard parse errors.
    Require a header-only final pathway report and an empty inference evidence
    list from this PGDB. Raw exports and archives remain unchanged.
    """
    flatpath = Path(flatpath)
    lines = [line.strip() for line in (flatpath / 'pathways.dat').read_text().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    if lines and lines != [EMPTY_EXPORT_ERROR]:
        return False  # Ordinary exports (including other errors) go to Camelot.
    reports = flatpath.parent / 'reports'
    summaries = sorted(reports.glob('pathways-report_*.txt'))
    if not summaries:
        raise ValueError('Cannot verify empty pathways.dat: missing final pathway report')
    # Multiple reports can reflect repeated inference. Fail closed if any differs.
    for report in summaries:
        rows = [line.strip() for line in report.read_text().splitlines()
                if line.strip() and not line.lstrip().startswith('#')]
        if len(rows) != 1 or tuple(part.strip() for part in rows[0].split('|')) != REPORT_COLUMNS:
            raise ValueError(f'Cannot verify empty pathways.dat: nonempty or invalid {report.name}')
    evidence = (reports / 'pwy-evidence-list.dat').read_text()
    if 'This file contains all inferred pathways and super-pathways' not in evidence or any(
        line.strip() and not line.lstrip().startswith(';;;')
        for line in evidence.splitlines()
    ):
        raise ValueError('Cannot verify empty pathways.dat: nonempty or invalid pathway evidence list')
    Path(outfile).write_text('\t'.join(PATHWAY_COLUMNS) + '\n')
    print('Pathway export: verified zero pathways; writing an empty pathway table', flush=True)
    return True
