"""Explicit sample/MAG scope for biological report tables."""
MAG_TABLES = {
    'contigs', 'orf_explorer', 'annotation_explorer', 'annotation_terms',
    'entities', 'contig_mags', 'mag_orf_explorer', 'mag_gene_explorer',
    'orf_groups', 'pathway_explorer', 'pathway_gene_explorer',
    'abundance_explorer', 'abundance', 'rna_explorer', 'rna_taxonomy_explorer',
    'mag_rna_explorer', 'orf_rna_explorer', 'rna_abundance_explorer', 'feature_explorer',
}
BUNDLE_TABLES = ('contigs', 'feature_explorer', 'orf_explorer', 'annotation_explorer',
                 'annotation_terms', 'rna_explorer', 'rna_taxonomy_explorer',
                 'abundance_explorer', 'pathway_explorer', 'pathway_gene_explorer',
                 'contig_mags', 'orf_groups')


def scope_clause(db, table, columns, scope):
    if not isinstance(scope, dict) or set(scope) - {'sample_id', 'entity_id'}:
        raise ValueError('Invalid sample/MAG scope')
    sample, mag = scope.get('sample_id'), scope.get('entity_id')
    if any(v is not None and (not isinstance(v, str) or not v or len(v) > 1000) for v in (sample, mag)):
        raise ValueError('Invalid sample/MAG identifier')
    if mag and not sample:
        raise ValueError('Choose a sample before choosing a MAG')
    if not sample:
        return '', []
    if not db.execute('SELECT 1 FROM samples WHERE sample_id=?', (sample,)).fetchone():
        raise ValueError('Unknown sample')
    if mag:
        entity = db.execute('SELECT entity_type FROM entities WHERE sample_id=? AND entity_id=?', (sample,mag)).fetchone()
        if not entity or entity[0] != 'MAG':
            raise ValueError('Unknown MAG in this sample')
        if not db.execute('SELECT 1 FROM contig_mags WHERE sample_id=? AND entity_id=?', (sample,mag)).fetchone():
            raise ValueError('Full contig-to-MAG membership is unavailable for this MAG')
        if table not in MAG_TABLES:
            raise ValueError('This table is run-wide metadata; choose a biological table for MAG filtering')
    if table == 'annotation_terms':
        membership = '' if not mag else ' AND EXISTS(SELECT 1 FROM contig_mags m WHERE m.sample_id=o.sample_id AND m.contig_id=o.contig_id AND m.entity_id=?)'
        return ('EXISTS(SELECT 1 FROM annotations a JOIN orfs o USING(sample_id,orf_id) WHERE a.annotation_id=t.annotation_id AND a.sample_id=?'+membership+')', [sample]+([mag] if mag else []))
    if 'sample_id' not in columns:
        raise ValueError('This table contains run-wide metadata and cannot be scoped by sample')
    clause, values = 't.sample_id=?', [sample]
    if not mag:
        return clause, values
    if 'entity_id' in columns:
        return clause+' AND t.entity_id=?', values+[mag]
    if table == 'abundance':
        expression = "m.contig_id=CASE WHEN t.feature_type='contig' THEN t.feature_id ELSE (SELECT o.contig_id FROM orfs o WHERE o.sample_id=t.sample_id AND o.orf_id=t.feature_id) END"
    elif table == 'orf_groups':
        expression = 'm.contig_id=(SELECT o.contig_id FROM orfs o WHERE o.sample_id=t.sample_id AND o.orf_id=t.member_orf_id)'
    elif 'contig_id' in columns:
        expression = 'm.contig_id=t.contig_id'
    else:
        raise ValueError('This table has no supported MAG membership relationship')
    return clause+' AND EXISTS(SELECT 1 FROM contig_mags m WHERE m.sample_id=t.sample_id AND m.entity_id=? AND '+expression+')', values+[mag]
