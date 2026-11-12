"""Scoped table and retained-file downloads, with an explicit availability manifest."""
import csv
import io
import json
import zipfile
from pathlib import Path
from urllib.parse import quote
from datetime import datetime, timezone
from metapathways.report_scope import MAG_TABLES, scope_clause


def selections(db, spec):
    if spec.get('mode') not in ('tables','complete'):
        raise ValueError('Choose tables or complete download mode')
    items = spec.get('selections')
    if spec.get('all_samples') is True:
        items = [{'sample_id':r[0]} for r in db.execute('SELECT sample_id FROM samples ORDER BY sample_id')]
    if not isinstance(items,list) or not items:
        raise ValueError('Choose at least one sample or MAG')
    unique = set()
    for scope in items:
        scope_clause(db,'contigs',['sample_id','contig_id'],scope)
        if not scope.get('sample_id'):
            raise ValueError('Each selection requires a sample')
        unique.add((scope['sample_id'],scope.get('entity_id')))
    # A whole sample already includes every MAG in it.
    whole = {s for s,m in unique if m is None}
    return [{'sample_id':s,**({'entity_id':m} if m else {})}
            for s,m in sorted(unique,key=lambda x:(x[0],x[1] or '')) if m is None or s not in whole]


def selected_files(db, scope):
    base = db.execute('SELECT output_path FROM samples WHERE sample_id=?',(scope['sample_id'],)).fetchone()[0]
    base_parts = Path(base).parts if base != '.' else ()
    mag = scope.get('entity_id')
    for row in db.execute('SELECT path,bytes,modified_utc FROM files ORDER BY path'):
        parts=Path(row[0]).parts
        if parts[:len(base_parts)] != base_parts: continue
        relative=parts[len(base_parts):]
        if not mag or relative[:3] == ('magsplitter','results',mag) or relative[:4] == ('results','pgdb','MAGs',mag):
            yield row


def checked_path(root, row):
    path = root/row[0]
    # A changed symlink must not turn a previously inventoried file into an escape.
    path.resolve().relative_to(root.resolve())
    if path.is_symlink() or any(p.is_symlink() for p in path.parents if p != root.parent):
        raise ValueError(f'Symbolic link in retained file path: {row[0]}')
    stat = path.stat()
    modified=datetime.fromtimestamp(stat.st_mtime,timezone.utc).isoformat()
    if stat.st_size != row[1] or modified != row[2]:
        raise ValueError(f'Output changed since report generation: {row[0]}; rebuild the report')
    return path


def preview(db, spec):
    scopes=selections(db,spec)
    return {'mode':spec['mode'],'selections':scopes,'items':[
        {'scope':scope,'retained_files':len(files),'retained_bytes':sum(r[1] for r in files)}
        for scope in scopes for files in [list(selected_files(db,scope))]],
        'note':'Complete MAG downloads include scoped sequence/GFF/GenBank subsets and MAG-specific files. Shared sample files are included only with whole samples. Missing files cannot be reconstructed.'}


def write_tables(db, archive, prefix, scope):
    from metapathways.reporting import VIEWS
    from metapathways.report_server import query, csv_value
    counts, sources = {}, set()
    for table in [*VIEWS,'sequence_explorer']:
        columns=[r[1] for r in db.execute('PRAGMA table_info("'+table+'")')]
        if scope.get('entity_id'):
            if table not in MAG_TABLES: continue
        elif 'sample_id' not in columns and table != 'annotation_terms': continue
        columns, cursor=query(db,{'table':table,'scope':scope},export=True,allow_sequence_scan=True)
        count=0
        with archive.open(prefix+'/tables/'+table+'.csv','w',force_zip64=True) as raw:
            with io.TextIOWrapper(raw,encoding='utf-8',newline='') as stream:
                writer=csv.writer(stream); writer.writerow(columns)
                for record in cursor:
                    writer.writerow([csv_value(v) for v in record]); count+=1
                    sources.update(record[i] for i,c in enumerate(columns) if c in ('source_id','taxonomy_source_id') and record[i] is not None)
        counts[table]=count
    with archive.open(prefix+'/tables/sources.csv','w',force_zip64=True) as raw:
        with io.TextIOWrapper(raw,encoding='utf-8',newline='') as stream:
            writer=csv.writer(stream); writer.writerow(['source_id','path','bytes','sha256','role'])
            for source in sorted(sources):
                writer.writerow([csv_value(v) for v in db.execute('SELECT * FROM sources WHERE source_id=?',(source,)).fetchone()])
    return counts


def write_mag_subsets(db, root, archive, prefix, scope):
    from metapathways.report_server import query
    sample, mag=scope['sample_id'],scope['entity_id']
    contigs={r[0] for r in db.execute('SELECT contig_id FROM contig_mags WHERE sample_id=? AND entity_id=?',(sample,mag))}
    aliases=set(contigs)
    alias_contigs={c:c for c in contigs}
    for c,original in db.execute('SELECT contig_id,original_id FROM contigs WHERE sample_id=?',(sample,)):
        if c in contigs and original and db.execute('SELECT COUNT(*) FROM contigs WHERE sample_id=? AND original_id=?',(sample,original)).fetchone()[0]==1:
            aliases.add(original); alias_contigs[original]=c
    counts={}
    seen_sequences={}
    # Stream original sequence records; record IDs in the table identify duplicates across sources.
    for kind in ('contig','protein','RNA','nucleotide'):
        _,records=query(db,{'table':'sequence_explorer','scope':scope,'filters':[{'column':'sequence_type','op':'eq','value':kind}]},export=True,allow_sequence_scan=True)
        count=0; seen=set()
        with archive.open(prefix+'/sequences/'+kind+'.fasta','w',force_zip64=True) as raw:
            with io.TextIOWrapper(raw,encoding='utf-8') as stream:
                for row in records:
                    stream.write('>'+str(row['sequence_record_id'])+'|'+row['description']+'\n')
                    seq=row['sequence']
                    for i in range(0,len(seq),80): stream.write(seq[i:i+80]+'\n')
                    count+=1; seen.add(row['contig_id'] if kind=='contig' else row['feature_id'])
        counts[kind+'_sequences']=count
        seen_sequences[kind]=seen
    base=db.execute('SELECT output_path FROM samples WHERE sample_id=?',(sample,)).fetchone()[0]
    baseparts=Path(base).parts if base != '.' else ()
    gbk_records=0; gff_records=0; gbk_contigs=set()
    for row in db.execute('SELECT path,bytes,modified_utc FROM files ORDER BY path'):
        parts=Path(row[0]).parts
        relative=parts[len(baseparts):]
        if parts[:len(baseparts)] != baseparts or not relative or relative[0]!='genbank': continue
        if not (row[0].endswith(('.gbk','.gb','.gbff','.annot.gff'))): continue
        path=checked_path(root,row)
        dest=prefix+'/annotations/'+quote(path.name,safe='')
        with archive.open(dest,'w',force_zip64=True) as raw, path.open() as source:
            with io.TextIOWrapper(raw,encoding='utf-8') as output:
                if path.suffix=='.gff':
                    output.write('##gff-version 3\n')
                    for line in source:
                        if line.startswith('##FASTA'): break
                        if not line.startswith('#') and line.split('\t')[0] in aliases:
                            output.write(line); gff_records+=1
                else:
                    keep=False; in_record=False
                    for line in source:
                        if line.startswith('LOCUS '):
                            if in_record: raise ValueError(f'Unterminated GenBank record in {path}')
                            in_record=True; record_id=line.split()[1]; keep=record_id in aliases
                        if keep: output.write(line)
                        if line.strip()=='//':
                            if keep: gbk_records+=1; gbk_contigs.add(alias_contigs[record_id])
                            in_record=False; keep=False
                    if in_record: raise ValueError(f'Unterminated GenBank record in {path}')
    counts.update(genbank_records=gbk_records,gff_features=gff_records)
    expected_genes={r[0] for r in db.execute('SELECT o.orf_id FROM orfs o JOIN contig_mags m USING(sample_id,contig_id) WHERE o.sample_id=? AND m.entity_id=? AND NOT EXISTS(SELECT 1 FROM rna_features r WHERE r.sample_id=o.sample_id AND r.rna_id=o.orf_id)',(sample,mag))}
    expected_rna={r[0] for r in db.execute('SELECT r.rna_id FROM rna_features r JOIN contig_mags m USING(sample_id,contig_id) WHERE r.sample_id=? AND m.entity_id=?',(sample,mag))}
    missing={'contig_sequences':contigs-seen_sequences['contig'], 'genbank_contigs':contigs-gbk_contigs,
             'protein_features':expected_genes-seen_sequences['protein'],
             'nucleotide_features':expected_genes-seen_sequences['nucleotide'], 'RNA_features':expected_rna-seen_sequences['RNA']}
    counts['missing_records']={key:{'count':len(ids),'examples':sorted(ids)[:100]} for key,ids in missing.items() if ids}
    return counts


def write_download(db, root, spec, destination):
    from metapathways.reporting import SCHEMA_VERSION
    scopes=selections(db,spec)
    manifest={'schema_version':SCHEMA_VERSION,'mode':spec['mode'],'selections':[]}
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for scope in scopes:
            prefix='sample-'+quote(scope['sample_id'],safe='')+'/'+('mag-'+quote(scope['entity_id'],safe='') if scope.get('entity_id') else 'complete-sample')
            entry={'scope':scope,'directory':prefix,'table_rows':write_tables(db,archive,prefix,scope),'files':[]}
            if spec['mode']=='complete':
                if scope.get('entity_id'):
                    entry['subset_records']=write_mag_subsets(db,root,archive,prefix,scope)
                    entry['availability_notes']=[f'No retained {k} available for this selection' for k,v in entry['subset_records'].items() if v==0]
                    entry['availability_notes'].extend(f'{data["count"]} missing {kind}; see subset_records.missing_records' for kind,data in entry['subset_records']['missing_records'].items())
                for row in selected_files(db,scope):
                    path=checked_path(root,row)
                    archive.write(path,prefix+'/files/'+row[0])
                    entry['files'].append({'path':row[0],'bytes':row[1]})
            manifest['selections'].append(entry)
        archive.writestr('manifest.json',json.dumps(manifest,indent=2))
        archive.writestr('README.txt',
            'MetaPathways selected results\n'
            'A complete sample includes all its indexed MAGs and retained inventoried files.\n'
            'Overlapping selections are deduplicated: a whole sample subsumes its MAGs.\n'
            'MAG membership uses the explicit sample/contig map, never protein taxonomy.\n'
            'MAG-only files include MAG-specific outputs and scoped FASTA/GFF/GenBank subsets.\n'
            'Shared BAMs, raw reads, logs and other sample-wide files require a whole-sample selection.\n'
            'Original input files outside the output directory and deleted intermediates are unavailable.\n'
            'Sequence CSVs contain full sequences; the single-record FASTA viewer uses bounded previews. FASTA headers prefix the sequence_record_id.\n'
            'GenBank records retain their original annotations; they are not reannotated.\n'
            'Empty tables/files mean no retained records, not biological absence. See manifest.json.\n'
            'Tables retain distinct row units and source IDs; do not sum duplicated measurements across views.\n')
