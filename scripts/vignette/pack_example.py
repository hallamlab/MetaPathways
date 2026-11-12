"""Validate and archive an already staged vitamin B5 example."""
import argparse,csv,hashlib,json,sqlite3,tarfile
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('example',type=Path);args=p.parse_args();root=args.example.resolve()
with sqlite3.connect(f'file:{root}/all/reports/results.sqlite?mode=ro',uri=True) as db:
    assert db.execute('pragma integrity_check').fetchone()[0]=='ok'
    assert db.execute('select count(*) from samples').fetchone()[0]==3
    assert db.execute("select count(*) from annotation_explorer where sample_id='Gastrointestinal_5' and orf_id='C2-G10'").fetchone()[0]==2
    assert db.execute("select covered_reactions,reactions from pathway_explorer where sample_id='Gastrointestinal_5' and entity_id='CAMI_OTU_97_68_0' and pathway_id='PANTO-PWY'").fetchone()==(2,4)
    assert db.execute("select count(*) from rna_taxonomy_explorer where sample_id='Gastrointestinal_5' and contig_id='Gastrointestinal_5-C2' and similarity=100").fetchone()[0]==2
    for path, in db.execute('select path from sources'):
        assert not Path(path).is_absolute() and (root/'all'/path).is_file(),path
for name in ['panc-annotations.csv','panto-pathway-genes.csv','panto-evidence.csv']:
    assert len(list(csv.DictReader((root/'exports'/name).open())))==2,name
files=[]
for path in sorted(root.rglob('*')):
    assert not path.is_symlink(),path
    if not path.is_file() or path==root/'manifest.json':continue
    rel=path.relative_to(root).as_posix()
    assert not any(part in ('MPDB','mp-containers','diagnostics','.metapathways') for part in Path(rel).parts),rel
    assert path.suffix not in ('.sif','.bam'),rel
    files.append(dict(path=rel,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
(root/'manifest.json').write_text(json.dumps(dict(snapshot=1,analysis_mp_version='4.0.8',report_mp_version='4.0.13',files=files),indent=2)+'\n')
archive=root.with_suffix('.tar.gz')
with tarfile.open(archive,'w:gz') as tar:tar.add(root,arcname=root.name)
digest=hashlib.sha256(archive.read_bytes()).hexdigest()
archive.with_name(archive.name+'.sha256').write_text(f'{digest}  {archive.name}\n')
print(f'Validated {len(files)} files; archive {archive.stat().st_size:,} bytes; sha256 {digest}')
