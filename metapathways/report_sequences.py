"""Stream retained FASTA/FASTQ sequences into the report, preserving source IDs."""
import gzip
import os
import re
from pathlib import Path

EXTENSIONS = {'.fa', '.fas', '.fsa', '.fasta', '.fna', '.faa', '.ffn', '.frn', '.fastq', '.fq'}


def records(path):
    opener = gzip.open if path.suffix.lower() == '.gz' else open
    with opener(path, 'rt') as stream:
        first = stream.readline()
        if not first:
            return
        if first.startswith('>'):
            header, parts = first[1:].strip(), []
            for line in stream:
                if line.startswith('>'):
                    yield header, ''.join(parts), None
                    header, parts = line[1:].strip(), []
                else:
                    parts.append(''.join(line.split()))
            yield header, ''.join(parts), None
        elif first.startswith('@'):
            while first:
                if not first.startswith('@'):
                    raise ValueError(f'Malformed FASTQ header in {path}')
                header, parts = first[1:].strip(), []
                for line in stream:
                    if line.startswith('+'):
                        break
                    parts.append(line.strip())
                else:
                    raise ValueError(f'Missing FASTQ quality separator in {path}')
                sequence = ''.join(parts)
                quality = ''
                while len(quality) < len(sequence):
                    line = stream.readline()
                    if not line:
                        raise ValueError(f'Truncated FASTQ quality in {path}')
                    quality += line.rstrip('\r\n')
                if len(quality) != len(sequence):
                    raise ValueError(f'FASTQ quality length differs from sequence in {path}')
                yield header, sequence, quality
                first = stream.readline()
        else:
            raise ValueError(f'Expected FASTA or FASTQ header in {path}')


def import_sequences(importer, sample, directory):
    db = importer.db
    for folder, dirs, files in os.walk(directory, followlinks=False):
        dirs[:] = sorted(d for d in dirs if not d.startswith('.') and d not in ('reports','work','conda-cache') and not (Path(folder)/d).is_symlink())
        for name in sorted(files):
            path = Path(folder)/name
            suffix = Path(path.stem).suffix.lower() if path.suffix.lower() == '.gz' else path.suffix.lower()
            if suffix not in EXTENSIONS or path.is_symlink():
                continue
            source_id = importer.source(path, 'retained sequence file')
            unmapped = 0
            for ordinal, (header, sequence, quality) in enumerate(records(path), 1):
                if not header or not sequence:
                    raise ValueError(f'Empty sequence identifier or sequence in {path}')
                identifier = header.split()[0]
                contig, feature, kind = None, None, 'read' if quality is not None else 'protein' if suffix == '.faa' else 'nucleotide'
                # Reads are not contigs, even if their identifiers happen to match.
                if quality is None:
                    if db.execute('SELECT 1 FROM contigs WHERE sample_id=? AND contig_id=?',(sample,identifier)).fetchone():
                        contig, kind = identifier, 'contig'
                    else:
                        originals=db.execute('SELECT contig_id FROM contigs WHERE sample_id=? AND original_id=?',(sample,identifier)).fetchall()
                        if len(originals)==1: contig, kind = originals[0][0], 'contig'
                    candidates = [identifier]
                    if '-C' in identifier:
                        candidates.append('C'+identifier.rsplit('-C', 1)[1])
                    matches=set()
                    for candidate in candidates:
                        for value in db.execute("SELECT contig_id,'CDS' FROM orfs WHERE sample_id=? AND orf_id=? AND contig_id IS NOT NULL UNION SELECT contig_id,feature_type FROM rna_features WHERE sample_id=? AND rna_id=?",(sample,candidate,sample,candidate)):
                            matches.add((candidate,*value))
                    # RNA abundance placeholders in the ORF table are not protein genes.
                    matches={m for m in matches if m[2]!='CDS' or not any(n[0]==m[0] and n[1]==m[1] and n[2]!='CDS' for n in matches)}
                    if contig is None and len(matches) == 1:
                        feature, contig, feature_kind = next(iter(matches))
                        if feature_kind in ('rRNA','tRNA','pseudogene'): kind = 'RNA'
                    coord = re.fullmatch(r'(.+):(\d+)-(\d+)\(([+-])\)', identifier.split('::')[-1])
                    if contig is None and coord:
                        c, start, end, strand = coord.groups()
                        match = db.execute('SELECT rna_id FROM rna_features WHERE sample_id=? AND contig_id=? AND start=? AND end=? AND strand=?', (sample,c,int(start)+1,int(end),strand)).fetchall()
                        if len(match) == 1: feature, contig, kind = match[0][0], c, 'RNA'
                if contig is None: unmapped += 1
                db.execute('INSERT INTO sequences(sample_id,sequence_id,record_number,sequence_type,contig_id,feature_id,description,length,sequence,quality,source_id) VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                           (sample,identifier,ordinal,kind,contig,feature,header,len(sequence),sequence,quality,source_id))
            if unmapped:
                importer.issue(sample,path.relative_to(importer.root),f'{unmapped} sequence records have no exact contig/feature link; available in sample exports, excluded from MAG-only exports.')
