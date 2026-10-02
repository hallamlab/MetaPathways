"""Attach real contig sequences to compact PathoLogic input in private staging."""
import csv
import json
import re
from pathlib import Path
from metapathways.pt_ec import normalize_ecs


def normalize_trna_name(name):
    """Separate MP's amino-acid label from a known anticodon; retain locus IDs."""
    return re.sub(
        r'(\.tRNA\d+-(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|'
        r'Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|fMet|SeC|Sec|Pyl))([ACGTU]{3})$',
        r'\1-\2', name)


def set_organism_taxon(inputs, taxon_id):
    """Apply an explicit NCBI taxon override to private PathoLogic input."""
    if not isinstance(taxon_id, int) or isinstance(taxon_id, bool) or taxon_id < 1:
        raise ValueError('NCBI taxon ID must be a positive integer')
    params = Path(inputs) / 'organism-params.dat'
    lines = params.read_text().splitlines()
    lines = [line for line in lines if line.split(None, 1)[:1] != ['NCBI-TAXON-ID']]
    lines.append(f'NCBI-TAXON-ID\t{taxon_id}')
    params.write_text('\n'.join(lines) + '\n')
    print(f'Pathway Tools organism taxon override: {taxon_id}', flush=True)


def prodigal_codon_tables(gff):
    """Read per-contig genetic codes from Prodigal's original model metadata."""
    tables, contig = {}, None
    if not gff.is_file():
        return tables
    for line in gff.read_text().splitlines():
        if line.startswith('# Sequence Data:'):
            match = re.search(r'seqhdr="([^" ]+)', line)
            contig = match.group(1) if match else None
        elif line.startswith('# Model Data:'):
            match = re.search(r'transl_table=(\d+)(?:;|$)', line)
            if contig is None or match is None:
                raise ValueError(f'Missing Prodigal contig/genetic-code metadata in {gff}')
            code = int(match.group(1))
            if code not in (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15):
                raise ValueError(f'Unsupported PathoLogic codon table {code} for {contig}')
            if contig in tables and tables[contig] != code:
                raise ValueError(f'Conflicting Prodigal codon tables for {contig}')
            tables[contig] = code
    return tables


def attach_sequences(inputs, sample_output):
    inputs, base = Path(inputs), Path(sample_output)
    sample = base.name
    table = base/'results/annotation_table'/f'{sample}.ptinput.tsv'
    fasta = base/'preprocessed'/f'{sample}.fasta'
    for file in (table, fasta):
        if not file.is_file():
            raise ValueError(f'Cannot prepare sequence-backed Pathway Tools input; missing {file}')
    with table.open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not {'orf_id', 'seqname', 'start', 'end', 'strand'} <= set(reader.fieldnames or []):
            raise ValueError(f'Pathway Tools annotation table lacks feature coordinates/contig/strand: {table}')
        features = {}
        for row in reader:
            if not row['orf_id'] or not row['seqname'] or row['orf_id'] in features:
                raise ValueError('Missing or duplicate feature/contig in annotation table: ' + row['orf_id'])
            start, end = int(row['start']), int(row['end'])
            if start < 1 or end < start or row['strand'] not in ('+', '-'):
                raise ValueError('Invalid source coordinates/strand for ' + row['orf_id'])
            if row['strand'] == '-':
                start, end = end, start
            features[row['orf_id']] = (row['seqname'], start, end)
    records, seen, corrected = {}, set(), 0
    normalized_ec_records = 0
    normalized_trna_names = {}
    texts, lines = [], []
    for line in (inputs/'0.pf').read_text().splitlines():
        if line.strip() == '//':
            texts.append('\n'.join(lines))
            lines = []
        else:
            lines.append(line)
    if any(line.strip() for line in lines):
        raise ValueError('Unterminated Pathway Tools feature record')
    for text in texts:
        if not text.strip():
            continue
        fields = dict(line.split('\t', 1) for line in text.strip().splitlines() if '\t' in line)
        identifier = fields.get('ID')
        if identifier in seen or identifier not in features:
            raise ValueError('Duplicate or unmapped Pathway Tools feature: ' + str(identifier))
        seen.add(identifier)
        contig, start, end = features[identifier]
        # MAG splitter expands representative annotations to member IDs but retains
        # the representative's coordinates. Restore each member's actual locus.
        corrected += (fields.get('STARTBASE'), fields.get('ENDBASE')) != (str(start), str(end))
        fixed = [line for line in text.strip().splitlines() if not line.startswith(('STARTBASE\t', 'ENDBASE\t'))]
        if fields.get('PRODUCT-TYPE', '').upper() == 'TRNA':
            old_name = fields.get('NAME', '')
            new_name = normalize_trna_name(old_name)
            if new_name != old_name:
                fixed = ['NAME\t' + new_name if line.startswith('NAME\t') else line for line in fixed]
                normalized_trna_names[identifier] = dict(original=old_name, staged=new_name)
        old_ecs = [line.split('\t', 1)[1] for line in fixed if line.startswith('EC\t')]
        ecs = normalize_ecs(old_ecs)
        normalized_ec_records += old_ecs != ecs
        fixed = [line for line in fixed if not line.startswith('EC\t')]
        fixed.extend('EC\t' + ec for ec in ecs)
        fixed.extend([f'STARTBASE\t{start}', f'ENDBASE\t{end}'])
        records.setdefault(contig, []).append(('\n'.join(fixed)+'\n//\n', start, end))
    if not records:
        raise ValueError('Pathway Tools input contains no feature records')
    sequences, name, parts = {}, None, []
    def save():
        if name in records:
            if name in sequences:
                raise ValueError('Duplicate contig sequence: ' + name)
            sequences[name] = ''.join(parts)
    with fasta.open() as stream:
        for line in stream:
            if line.startswith('>'):
                save()
                name, parts = line[1:].split()[0], []
            elif name in records:
                parts.append(line.strip())
    save()
    missing = set(records)-set(sequences)
    if missing:
        raise ValueError('Missing Pathway Tools contig sequences: ' + ', '.join(sorted(missing)[:10]))
    # Validate everything before replacing the genetic-element manifest.
    for contig, features in records.items():
        for _, start, end in features:
            if min(start, end) < 1 or max(start, end) > len(sequences[contig]):
                raise ValueError(f'Pathway Tools coordinates {start}..{end} outside {contig} (length {len(sequences[contig])})')
    codon_tables = prodigal_codon_tables(base/'orf_prediction'/f'{sample}.cds.gff')
    elements = []
    for number, contig in enumerate(sorted(records), 1):
        identifier = f'contig_{number}'
        (inputs/(identifier+'.pf')).write_text(''.join(record[0] for record in records[contig]))
        seq = sequences[contig]
        (inputs/(identifier+'.fasta')).write_text('>'+identifier+'\n' + '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80))+'\n')
        code = f'CODON-TABLE\t{codon_tables[contig]}\n' if contig in codon_tables else ''
        elements.append(f'ID\t{identifier}\nNAME\t{contig}\nTYPE\t:CONTIG\n{code}ANNOT-FILE\t{identifier}.pf\nSEQ-FILE\t{identifier}.fasta\n//\n')
    (inputs/'genetic-elements.dat').write_text(''.join(elements))
    (inputs/'sequence-input.json').write_text(json.dumps(dict(contigs=len(sequences), features=len(seen),
        restored_coordinates=corrected, normalized_ec_records=normalized_ec_records,
        normalized_trna_names=normalized_trna_names,
        codon_tables={c: codon_tables[c] for c in records if c in codon_tables}, fasta=str(fasta), feature_table=str(table)), indent=2)+'\n')
    if normalized_trna_names:
        print(f'Separated anticodon triplets in {len(normalized_trna_names)} tRNA names; feature IDs and original inputs retained', flush=True)
    if normalized_ec_records:
        print(f'Normalized EC entries in {normalized_ec_records} features; original inputs retained', flush=True)
    print(f'Attached {len(sequences)} contig sequences to {len(seen)} Pathway Tools features; restored {corrected} feature coordinates', flush=True)
