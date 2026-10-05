"""Readable aliases for explicit Pathway Tools organism taxon overrides."""
import argparse

SCOPES = {'all': 131567, 'bacteria': 2, 'archaea': 2157, 'eukaryotes': 2759}
ALIASES = {'euks': 'eukaryotes'}


def positive_taxon(value):
    try:
        taxon = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError('taxon ID must be a positive integer')
    if taxon < 1:
        raise argparse.ArgumentTypeError('taxon ID must be a positive integer')
    return taxon


def scope_name(value):
    if value == 'prokaryotes':
        raise argparse.ArgumentTypeError(
            'prokaryotes combines Bacteria and Archaea and has no supported single '
            'NCBI taxon here; choose bacteria, archaea, or all (includes eukaryotes)')
    return ALIASES.get(value, value)


def add_taxonomy_options(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--taxon_id', type=positive_taxon,
                       help='Override PGDB NCBI taxon in private inputs; applies to every selected entity')
    group.add_argument('--taxonomic_scope', type=scope_name, choices=tuple(SCOPES),
                       help='Named PGDB taxon override; all means cellular life; euks aliases eukaryotes. '
                            'Taxonomic pruning is enabled by default. '
                            'Default: all (cellular life).')


def resolve_taxon(args):
    scope = getattr(args, 'taxonomic_scope', None)
    taxon = getattr(args, 'taxon_id', None)
    if scope and taxon is not None:
        raise ValueError('--taxon_id and --taxonomic_scope are mutually exclusive')
    return SCOPES[scope] if scope else (taxon if taxon is not None else SCOPES['all'])


def add_pruning_options(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--taxprune', dest='taxprune', action='store_true', default=True,
                       help='Enable taxonomic pruning [default]')
    group.add_argument('--no_taxprune', dest='taxprune', action='store_false',
                       help='Disable taxonomic pruning and perform unpruned rescoring')
