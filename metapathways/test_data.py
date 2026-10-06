"""Copy the bundled test inputs and reference seeds into a writable workspace."""
import argparse
from pathlib import Path
import shutil

SEEDS = ('functional/swissprot_test', 'taxonomic/SILVA_SSU_test',
         'taxonomic/SILVA_LSU_test', 'test_reference.json')


def prepare(output, fixtures=None):
    fixtures = Path(fixtures) if fixtures else Path(__file__).parent / 'regtests'
    root = Path(output).expanduser().resolve()
    inputs = fixtures / 'cami_test'
    files = [(p, root / 'cami-test' / p.relative_to(inputs))
             for p in sorted(inputs.rglob('*')) if p.is_file()]
    if not files or not (inputs / 'all.tsv').is_file():
        raise ValueError('The installed test bundle is incomplete; reinstall MetaPathways.')
    files.extend((fixtures / 'test_db' / name, root / 'MPDB' / name) for name in SEEDS)
    # Preflight the complete copy before writing anything. Never follow a link into
    # an installation or replace a user's modified inputs/reference files.
    for source, target in files:
        if not source.is_file():
            raise ValueError(f'Missing bundled test file: {source}')
        for p in (target, *target.parents):
            if p == root:
                break
            if p.is_symlink():
                raise ValueError(f'Test destination contains a symbolic link: {p}; choose a new directory.')
        if any(p.exists() and not p.is_dir() for p in target.parents):
            raise ValueError(f'Test destination parent is not a directory: {target.parent}')
        if target.exists() and (not target.is_file() or target.read_bytes() != source.read_bytes()):
            raise ValueError(f'Test file differs: {target}; choose a new directory.')
    for source, target in files:
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            # Exclusive creation also protects against a concurrent preparation.
            with source.open('rb') as src, target.open('xb') as dst:
                shutil.copyfileobj(src, dst)
    return root


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--output_dir', required=True,
                        help='workspace for cami-test/ inputs and MPDB/ reference seeds')
    args = parser.parse_args(argv)
    root = prepare(args.output_dir)
    print(f'Test workspace: {root}')
    print('Inputs: cami-test/all.tsv (three samples), single.tsv, pair.tsv')
    print('Next: change to this workspace and run metapathways build_db --test -d MPDB')


if __name__ == '__main__':
    main()
