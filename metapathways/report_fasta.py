"""Exact, one-record FASTA inspection; bulk sequences belong in downloads."""
LIMIT = 10000


def options(db, spec):
    if not isinstance(spec,dict) or set(spec)-{'sample_id','feature_id','contig_id'}:
        raise ValueError('Invalid FASTA lookup')
    sample=spec.get('sample_id')
    targets=[k for k in ('feature_id','contig_id') if spec.get(k)]
    if not isinstance(sample,str) or not sample or len(targets)!=1:
        raise ValueError('Choose one feature or contig in one sample')
    key=targets[0]; value=spec[key]
    if not isinstance(value,str) or not value or len(value)>1000:
        raise ValueError('Invalid feature/contig identifier')
    clause='q.'+key+'=?'
    if key=='contig_id': clause+=" AND q.sequence_type='contig'"
    rows=db.execute('SELECT q.sequence_record_id,q.sequence_id,q.sequence_type,q.description,q.length,s.path AS source FROM sequences q JOIN sources s USING(source_id) WHERE q.sample_id=? AND '+clause+' ORDER BY q.sequence_type,s.path,q.record_number LIMIT 501',(sample,value)).fetchall()
    return {'records':[dict(r) for r in rows[:500]],'more_records':len(rows)>500,
            'note':'Choose one retained sequence. No sequence is translated or inferred here. Bulk records are available from Downloads.'}


def record(db, spec, full=False):
    if not isinstance(spec,dict) or set(spec)-{'sample_id','sequence_record_id'}:
        raise ValueError('Invalid FASTA record request')
    sample,identifier=spec.get('sample_id'),spec.get('sequence_record_id')
    if not isinstance(sample,str) or not sample or not isinstance(identifier,int) or isinstance(identifier,bool) or identifier<1:
        raise ValueError('Choose one sequence record in one sample')
    sequence_column='sequence' if full else f'substr(sequence,1,{LIMIT})'
    row=db.execute('SELECT description,length,'+sequence_column+' AS sequence FROM sequences WHERE sample_id=? AND sequence_record_id=?',(sample,identifier)).fetchone()
    if row is None: raise ValueError('Sequence record is unavailable in this sample')
    sequence=row['sequence']
    fasta='>'+row['description']+'\n'+'\n'.join(sequence[i:i+80] for i in range(0,len(sequence),80))+'\n'
    return {'fasta':fasta,'length':row['length'],'truncated':len(sequence)<row['length']}
