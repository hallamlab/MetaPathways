"""Taxonomy provenance for individual protein hits and within-database LCA."""
import re

NOT_COMPUTED = 'Not computed'
UNCLASSIFIED = 'Unclassified'


def supports_taxonomy(database):
    return any(name in database.lower() for name in ('swissprot', 'uniref', 'eggnog'))


def hit_taxid(hit, database):
    database = database.lower()
    if 'eggnog' in database:
        match = re.match(r'^(\d+)\.', hit.get('target', ''))
    else:
        field = 'OX' if 'swissprot' in database else 'TaxID' if 'uniref' in database else None
        if field is None:
            return None
        # Parsed reports historically replace '=' with spaces; accept both forms.
        text = ' '.join(str(hit.get(k, '')) for k in ('product', 'comment'))
        match = re.search(r'\b' + field + r'(?:\s*=\s*|\s+)(\d+)(?=\s|$)', text)
    return match.group(1) if match else None


def taxon_label(lca, taxid):
    if taxid is None or str(taxid) not in lca.taxid_to_ptaxid:
        return UNCLASSIFIED
    taxid = str(taxid)
    return lca.get_preferred_taxonomy(taxid) or lca.id_to_name[taxid]


def hit_taxonomy(hit, database, lca):
    if not supports_taxonomy(database):
        return '', NOT_COMPUTED
    taxid = hit_taxid(hit, database)
    return taxid or '', taxon_label(lca, taxid)


def raw_lca(hits, database, lca):
    """Only valid IDs from score-qualified hits within this database enter LCA."""
    eligible = [h for h in hits if h['bitscore'] >= lca.lca_min_score]
    if not eligible:
        return None
    threshold = max(h['bitscore'] for h in eligible) * (1 - lca.lca_top_percent / 100)
    ids = {hit_taxid(h, database) for h in eligible if h['bitscore'] >= threshold}
    ids = sorted(i for i in ids if i in lca.taxid_to_ptaxid)
    if not ids:
        return None
    return str(lca.getTaxonomy(ids, taxid=True, return_id=True))
