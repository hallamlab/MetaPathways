"""Import RNA loci and their separate reference-hit evidence for reporting."""
import csv
import re
from urllib.parse import unquote


def import_rna(importer, sample, directory):
    db = importer.db
    # The final GFF uses the identifiers used by abundance and downstream outputs.
    annotated = list((directory/'genbank').glob('*.annot.gff'))
    paths = annotated or (list((directory/'orf_prediction').glob('*.rRNA.gff')) +
                          list((directory/'orf_prediction').glob('*.tRNA.gff')))
    if len(annotated) > 1:
        raise ValueError(f'Multiple annotated GFF files in {directory}/genbank')
    for path in sorted(paths):
        source = importer.source(path, 'RNA feature coordinates and annotations')
        with path.open() as stream:
            for line in stream:
                if line.startswith('##FASTA'):
                    break
                if not line.strip() or line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) != 9:
                    raise ValueError(f'Malformed GFF row in {path}')
                attrs = dict((unquote(k), unquote(v)) for k, v in
                             (item.split('=', 1) for item in f[8].split(';') if '=' in item))
                kind = f[2]
                if kind not in ('rRNA', 'tRNA') and not (kind == 'pseudogene' and 'anticodon' in attrs):
                    continue
                contig, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
                name = attrs.get('product') or attrs.get('Name') or ''
                subtype = re.search(r'(?:5\.8|5|16|18|23|28)S', name)
                identifier = attrs.get('ID') or f'{contig}:{start}-{end}({strand}):{kind}'
                db.execute('INSERT OR IGNORE INTO contigs(sample_id,contig_id) VALUES(?,?)', (sample, contig))
                db.execute('INSERT INTO rna_features VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)',
                           (sample, identifier, contig, kind, subtype.group() if subtype else None,
                            start, end, strand, name or None, attrs.get('isotype') or attrs.get('isotype_prediction'),
                            attrs.get('anticodon'), f[5] if f[5] != '.' else None, source))
    for path in sorted((directory/'results/rRNA').glob('*.rRNA.stats.txt')):
        source = importer.source(path, 'rRNA reference taxonomy')
        reference = path.name.removesuffix('.rRNA.stats.txt')
        prefix = directory.name + '.'
        if reference.startswith(prefix):
            reference = reference[len(prefix):]
        with path.open() as stream:
            reader = csv.reader(stream, delimiter='\t')
            for header in reader:
                if header[:3] == ['sequence', 'start', 'end']:
                    break
            else:
                importer.issue(sample, path.relative_to(importer.root), 'No rRNA taxonomy header found.')
                continue
            for row in reader:
                if not row:
                    continue
                if len(row) != 7:
                    raise ValueError(f'Expected seven rRNA taxonomy fields in {path}')
                query, hit_start, hit_end, similarity, evalue, bitscore, taxonomy = row
                if hit_start == '-':
                    continue
                # bedtools coordinates are zero-based; hit offsets are within the
                # extracted RNA, never coordinates on the assembled contig.
                match = re.search(r'^(.+):(\d+)-(\d+)\(([+-])\)$', query.split('::')[-1])
                matches = []
                if match:
                    contig, start, end, strand = match.groups()
                    matches = db.execute('SELECT rna_id FROM rna_features WHERE sample_id=? AND contig_id=? AND start=? AND end=? AND strand=? AND feature_type=?',
                                         (sample, contig, int(start)+1, int(end), strand, 'rRNA')).fetchall()
                else:
                    matches = db.execute('SELECT rna_id FROM rna_features WHERE sample_id=? AND rna_id=? AND feature_type=?', (sample, query, 'rRNA')).fetchall()
                identifier = matches[0][0] if len(matches) == 1 else None
                if identifier is None:
                    importer.issue(sample, path.relative_to(importer.root), f'RNA taxonomy query {query}: no unique RNA locus match; evidence retained without a locus link.')
                db.execute('INSERT INTO rna_taxonomy(sample_id,rna_id,query_id,reference_db,hit_start,hit_end,similarity,evalue,bitscore,taxonomy,source_id) VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                           (sample, identifier, query, reference, int(hit_start), int(hit_end), float(similarity), float(evalue), float(bitscore), taxonomy, source))
