#!/usr/bin/env python3
"""Prepare tiny real CAMI inputs; requires minimap2, not MetaPathways or Pathway Tools.

Select three abundant genomes per sample, retain up to 50 kb of one contig each,
then retain intact original read pairs with a primary MAPQ >=20 alignment to
those regions from the first 250,000 pairs. This is deliberately coverage-biased
interface test data, never an abundance or accuracy benchmark.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

SAMPLES = ('Urogenital_22', 'Gastrointestinal_5', 'Skin_28')
FIELDS = ('sample_id', 'assembly', 'read_layout', 'reads_1', 'reads_2', 'mag_map')


def fasta(path):
    with open(path) as handle:
        name, seq = None, []
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, ''.join(seq)


def pair_records(handle):
    while True:
        a = [handle.readline() for _ in range(4)]
        if not a[0]:
            return
        b = [handle.readline() for _ in range(4)]
        for record in (a, b):
            if (not record[0].startswith('@') or not record[2].startswith('+')
                    or len(record[1].strip()) != len(record[3].strip()) or not record[3]):
                raise ValueError('Invalid/truncated FASTQ')
        if not a[0].split()[0].endswith('/1') or a[0].split()[0][:-2] != b[0].split()[0][:-2] or not b[0].split()[0].endswith('/2'):
            raise ValueError('Source FASTQ is not adjacent /1, /2 pairs')
        yield a, b


def write_gzip(path, text):
    # No timestamp or source filename: reproducible compressed payload.
    with open(path, 'wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as handle:
            handle.write(text.encode())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-manifest', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--minimap2', default='minimap2')
    args = parser.parse_args()
    out = Path(args.output)
    if out.exists():
        parser.error('Output already exists; choose a new directory')
    out.mkdir(parents=True)
    for part in ('assemblies', 'reads', 'mag_maps'):
        (out / 'inputs' / part).mkdir(parents=True)
    with open(args.source_manifest) as handle:
        sources = {r['sample_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    rows, stats = [], []
    with tempfile.TemporaryDirectory(prefix='mp-cami-test-') as tmp:
        tmp = Path(tmp)
        for sample in SAMPLES:
            source = sources[sample]
            if source['read_layout'] != 'interleaved' or source.get('reads_2'):
                raise ValueError(f'{sample}: expected original interleaved CAMI reads')
            assembly = Path(source['assembly'])
            mapping = assembly.parent / 'gsa_mapping.tsv'
            with mapping.open() as handle:
                candidates = [r for r in csv.DictReader(handle, delimiter='\t')
                              if int(r['end_position']) - int(r['start_position']) + 1 >= 20000]
            candidates.sort(key=lambda r: (-int(r['number_reads']) / (int(r['end_position']) - int(r['start_position']) + 1), r['#anonymous_contig_id']))
            chosen, genomes = {}, set()
            for r in candidates:
                if r['genome_id'] not in genomes:
                    chosen[r['#anonymous_contig_id']] = r
                    genomes.add(r['genome_id'])
                if len(chosen) == 3:
                    break
            if len(chosen) != 3:
                raise ValueError(f'{sample}: fewer than three suitable genomes')
            sequences = {}
            for name, seq in fasta(assembly):
                if name in chosen:
                    sequences[name] = seq[:50000]
                if len(sequences) == len(chosen):
                    break
            if sequences.keys() != chosen.keys():
                raise ValueError(f'{sample}: missing selected contigs')
            fa = ''.join(f'>{name}\n{sequences[name]}\n' for name in sorted(sequences))
            (tmp / 'assembly.fasta').write_text(fa)
            with gzip.open(source['reads_1'], 'rt') as inp, (tmp / 'r1.fq').open('w') as r1, (tmp / 'r2.fq').open('w') as r2:
                scanned = 0
                for a, b in pair_records(inp):
                    r1.writelines(a)
                    r2.writelines(b)
                    scanned += 1
                    if scanned == 250000:
                        break
            cmd = [args.minimap2, '-ax', 'sr', '-t', '2', str(tmp / 'assembly.fasta'), str(tmp / 'r1.fq'), str(tmp / 'r2.fq')]
            with (tmp / 'align.sam').open('w') as sam, (out / f'{sample}.selection.log').open('w') as log:
                subprocess.run(cmd, stdout=sam, stderr=log, check=True)
            hits, per_contig = {}, {name: set() for name in sequences}
            with (tmp / 'align.sam').open() as sam:
                for line in sam:
                    if line.startswith('@'):
                        continue
                    f = line.split('\t')
                    flag = int(f[1])
                    if flag & (4 | 256 | 2048) or int(f[4]) < 20:
                        continue
                    name = re.sub(r'/[12]$', '', f[0])
                    hits.setdefault(name, set()).add(f[2])
            r1_out, r2_out, kept = [], [], 0
            with (tmp / 'r1.fq').open() as r1, (tmp / 'r2.fq').open() as r2:
                for _ in range(scanned):
                    a, b = [r1.readline() for _ in range(4)], [r2.readline() for _ in range(4)]
                    name = a[0].split()[0][1:-2]
                    if name in hits:
                        r1_out.extend(a)
                        r2_out.extend(b)
                        kept += 1
                        for contig in hits[name]:
                            per_contig[contig].add(name)
                        if kept == 3000:
                            break
            if any(len(v) < 10 for v in per_contig.values()):
                raise ValueError(f'{sample}: insufficient matched pairs: {[(k, len(v)) for k, v in per_contig.items()]}')
            write_gzip(out / 'inputs' / 'assemblies' / f'{sample}.fasta.gz', fa)
            write_gzip(out / 'inputs' / 'reads' / f'{sample}_R1.fastq.gz', ''.join(r1_out))
            write_gzip(out / 'inputs' / 'reads' / f'{sample}_R2.fastq.gz', ''.join(r2_out))
            (out / 'inputs' / 'mag_maps' / f'{sample}.tsv').write_text(''.join(f'{name}\tCAMI_{re.sub("[^A-Za-z0-9_]", "_", chosen[name]["genome_id"])}\n' for name in sorted(sequences)))
            rows.append(dict(zip(FIELDS, (sample, f'inputs/assemblies/{sample}.fasta.gz', 'paired', f'inputs/reads/{sample}_R1.fastq.gz', f'inputs/reads/{sample}_R2.fastq.gz', f'inputs/mag_maps/{sample}.tsv'))))
            stats.append({'sample_id': sample, 'source_assembly': str(assembly), 'source_reads': source['reads_1'], 'source_mapping': str(mapping), 'scanned_pairs': scanned, 'retained_pairs': kept, 'assembly_bases': sum(map(len, sequences.values())), 'regions': [{**chosen[name], 'retained_contig_start_1based': 1, 'retained_contig_end_1based': len(sequences[name]), 'aligned_pairs': len(per_contig[name])} for name in sorted(sequences)]})
            print(f'{sample}: {stats[-1]["assembly_bases"]} bases, {kept} pairs, 3 genome bins', flush=True)
    for label, subset in [('single', rows[:1]), ('pair', rows[1:]), ('all', rows)]:
        with (out / f'{label}.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(subset)
    provenance = {'source_dataset': 'CAMI II human-associated short-read gold-standard assemblies and simulated reads', 'source_doi': '10.4126/FRL01-006425518', 'selection': __doc__, 'minimap2_version': subprocess.check_output([args.minimap2, '--version'], text=True).strip(), 'samples': stats, 'sha256': {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}}
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
