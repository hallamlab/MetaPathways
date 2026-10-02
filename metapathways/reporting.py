"""Read existing MP results into a sample-scoped, queryable report database."""
import csv
from contextlib import closing
import fcntl
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile
from datetime import datetime, timezone
from urllib.parse import quote

SCHEMA_VERSION = 1
SCHEMA = '''
PRAGMA foreign_keys=ON;
CREATE TABLE samples(sample_id TEXT PRIMARY KEY, output_path TEXT NOT NULL);
CREATE TABLE sources(source_id INTEGER PRIMARY KEY, path TEXT UNIQUE, bytes INTEGER, sha256 TEXT, role TEXT);
CREATE TABLE contigs(sample_id TEXT, contig_id TEXT, original_id TEXT, length INTEGER,
 PRIMARY KEY(sample_id,contig_id), FOREIGN KEY(sample_id) REFERENCES samples);
CREATE TABLE orfs(sample_id TEXT, orf_id TEXT, contig_id TEXT, length INTEGER, start INTEGER, end INTEGER,
 strand TEXT, target TEXT, product TEXT, taxonomy TEXT, annotation_present INTEGER NOT NULL DEFAULT 0,
 PRIMARY KEY(sample_id,orf_id), FOREIGN KEY(sample_id,contig_id) REFERENCES contigs);
CREATE TABLE annotations(annotation_id INTEGER PRIMARY KEY, sample_id TEXT, orf_id TEXT,
 reference_db TEXT, target TEXT, product TEXT, score REAL, ec TEXT, reaction TEXT, source_id INTEGER,
 FOREIGN KEY(sample_id,orf_id) REFERENCES orfs, FOREIGN KEY(source_id) REFERENCES sources);
CREATE TABLE annotation_terms(annotation_id INTEGER, term_type TEXT, term TEXT,
 PRIMARY KEY(annotation_id,term_type,term), FOREIGN KEY(annotation_id) REFERENCES annotations);
CREATE TABLE entities(sample_id TEXT, entity_id TEXT, entity_type TEXT, pathway_status TEXT, last_task_status TEXT,
 PRIMARY KEY(sample_id,entity_id), FOREIGN KEY(sample_id) REFERENCES samples);
CREATE TABLE contig_mags(sample_id TEXT, contig_id TEXT, entity_id TEXT, original_mag_id TEXT, source_id INTEGER,
 PRIMARY KEY(sample_id,contig_id,entity_id), FOREIGN KEY(sample_id,contig_id) REFERENCES contigs,
 FOREIGN KEY(sample_id,entity_id) REFERENCES entities, FOREIGN KEY(source_id) REFERENCES sources);
CREATE INDEX original_contigs ON contigs(sample_id,original_id);
CREATE TABLE entity_orfs(sample_id TEXT, entity_id TEXT, orf_id TEXT,
 PRIMARY KEY(sample_id,entity_id,orf_id), FOREIGN KEY(sample_id,entity_id) REFERENCES entities,
 FOREIGN KEY(sample_id,orf_id) REFERENCES orfs);
CREATE TABLE orf_groups(sample_id TEXT, representative_orf_id TEXT, member_orf_id TEXT,
 PRIMARY KEY(sample_id,representative_orf_id,member_orf_id),
 FOREIGN KEY(sample_id,representative_orf_id) REFERENCES orfs,
 FOREIGN KEY(sample_id,member_orf_id) REFERENCES orfs);
CREATE TABLE pathways(sample_id TEXT, entity_id TEXT, pathway_id TEXT, name TEXT, score REAL,
 reactions INTEGER, covered_reactions INTEGER, reported_orf_count INTEGER, source_id INTEGER,
 PRIMARY KEY(sample_id,entity_id,pathway_id), FOREIGN KEY(sample_id,entity_id) REFERENCES entities,
 FOREIGN KEY(source_id) REFERENCES sources);
CREATE TABLE pathway_orfs(sample_id TEXT, entity_id TEXT, pathway_id TEXT, orf_id TEXT,
 PRIMARY KEY(sample_id,entity_id,pathway_id,orf_id),
 FOREIGN KEY(sample_id,entity_id,pathway_id) REFERENCES pathways,
 FOREIGN KEY(sample_id,orf_id) REFERENCES orfs);
CREATE TABLE abundance(sample_id TEXT, feature_type TEXT, feature_id TEXT, measurement TEXT, value REAL,
 source_id INTEGER, PRIMARY KEY(sample_id,feature_type,feature_id,measurement),
 FOREIGN KEY(sample_id) REFERENCES samples, FOREIGN KEY(source_id) REFERENCES sources);
CREATE TABLE issues(issue_id INTEGER PRIMARY KEY, sample_id TEXT, source TEXT, message TEXT);
CREATE TABLE files(path TEXT PRIMARY KEY, bytes INTEGER, modified_utc TEXT);
CREATE TABLE execution(run_id TEXT, command TEXT, task TEXT, label TEXT, status TEXT,
 elapsed_seconds REAL, error TEXT, source TEXT);
CREATE INDEX annotations_orf ON annotations(sample_id,orf_id);
CREATE INDEX annotations_reference ON annotations(reference_db);
CREATE INDEX terms_value ON annotation_terms(term_type,term);
CREATE INDEX entity_orfs_orf ON entity_orfs(sample_id,orf_id);
CREATE INDEX pathway_orfs_orf ON pathway_orfs(sample_id,orf_id);
CREATE INDEX orfs_contig ON orfs(sample_id,contig_id);
CREATE VIEW orf_explorer AS
 SELECT o.*, c.original_id AS original_contig_id, c.length AS contig_length
 FROM orfs o LEFT JOIN contigs c USING(sample_id,contig_id);
CREATE VIEW annotation_explorer AS
 SELECT a.*, o.contig_id, o.taxonomy, c.original_id AS original_contig_id
 FROM annotations a JOIN orfs o USING(sample_id,orf_id)
 LEFT JOIN contigs c USING(sample_id,contig_id);
CREATE VIEW pathway_explorer AS
 SELECT p.*, e.entity_type,
 (SELECT COUNT(*) FROM pathway_orfs g WHERE g.sample_id=p.sample_id AND g.entity_id=p.entity_id
 AND g.pathway_id=p.pathway_id) AS linked_orf_count
 FROM pathways p JOIN entities e USING(sample_id,entity_id);
CREATE VIEW pathway_gene_explorer AS
 SELECT g.*, p.name AS pathway_name, p.score AS pathway_score, e.entity_type,
 o.contig_id, o.product, o.taxonomy, o.annotation_present
 FROM pathway_orfs g JOIN pathways p USING(sample_id,entity_id,pathway_id)
 JOIN orfs o USING(sample_id,orf_id) JOIN entities e USING(sample_id,entity_id);
CREATE VIEW mag_orf_explorer AS
 SELECT m.sample_id,m.entity_id,m.contig_id,o.orf_id,o.product,o.taxonomy,o.annotation_present,m.source_id
 FROM contig_mags m JOIN orfs o USING(sample_id,contig_id);
CREATE VIEW mag_gene_explorer AS
 SELECT m.*, e.entity_type, o.contig_id, o.product, o.taxonomy, o.annotation_present
 FROM entity_orfs m JOIN entities e USING(sample_id,entity_id) JOIN orfs o USING(sample_id,orf_id);
'''
VIEWS = {
 'samples': ('Samples', 'One row per sample output directory.'),
 'contigs': ('Contigs', 'One row per sample and contig; original identifiers come from the mapping file.'),
 'orf_explorer': ('ORFs and taxonomy', 'One row per sample and ORF; the primary annotation and reported taxonomy are preserved.'),
 'annotation_explorer': ('Functional annotations', 'One row per reference annotation record; an ORF can have multiple records.'),
 'annotation_terms': ('EC and reaction terms', 'One row per annotation and EC/reaction term; join by annotation_id.'),
 'entities': ('Communities and MAGs', 'Pathway output availability and recorded task status; unavailable is not biological absence.'),
 'contig_mags': ('Contig-to-MAG membership', 'Explicit original contig-to-MAG assignments, joined through the contig identifier map.'),
 'mag_orf_explorer': ('MAG ORFs', 'All reported ORFs on explicitly mapped MAG contigs; requires the preserved contig-to-MAG map.'),
 'mag_gene_explorer': ('MAG input genes', 'Only genes explicitly listed in each MAG Pathway Tools input; not complete MAG membership.'),
 'orf_groups': ('Collapsed ORF groups', 'Explicit representative/member mappings from ptools/orf_map.txt; no membership is inferred.'),
 'pathway_explorer': ('Pathways', 'One row per sample, community/MAG and pathway. Reported scores are not recalculated.'),
 'pathway_gene_explorer': ('Pathway genes', 'One row per explicit pathway/ORF association, without multiplying by reference hits.'),
 'abundance': ('Read abundance', 'One row per feature and original measurement column. Values are copied, not recalculated or validated.'),
 'execution': ('Execution history', 'Every retained invocation; repeated or reused tasks are not independent biological results.'),
 'issues': ('Import notes', 'Missing files, unmatched identifiers and other limitations detected while reading outputs.'),
 'sources': ('Indexed sources', 'Relative paths and SHA-256 checksums of the files used to build the database.'),
 'files': ('All output files', 'Output inventory, excluding report products and temporary workflow state.'),
}


def number(value, integer=False):
    if value in (None, '', 'nan', 'NA', 'None'):
        return None
    result = float(value)
    if not __import__('math').isfinite(result):
        return None
    return int(result) if integer else result


def rows(path):
    """TSV header normalization; preserve blanks and literal identifier strings."""
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if reader.fieldnames:
            reader.fieldnames = [h.lstrip('# ').strip() for h in reader.fieldnames]
        for line, row in enumerate(reader, 2):
            if None in row:
                raise ValueError(f'{path}:{line}: more values than header columns')
            yield {k: v if v is not None else '' for k, v in row.items()}


def sample_dirs(root):
    def is_sample(p):
        return any((p / name).is_dir() for name in ('results/annotation_table', 'preprocessed', 'ptools'))
    if is_sample(root):
        return [root]
    return sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith('.') and is_sample(p))


class Importer:
    def __init__(self, db, root):
        self.db, self.root = db, root
        self.source_ids = {}
        self.source_stats = {}
        self.pgdb_status = {}
        # Ignore incomplete or stale pathway tables from a failed/skipped retry.
        for summary in sorted(root.glob('logs/analysis_wf/*/summary.json'), key=lambda p: p.stat().st_mtime_ns):
            for task in json.loads(summary.read_text()).get('tasks', []):
                if task.get('sample') and task.get('entity'):
                    self.pgdb_status[(task['sample'], task['entity'])] = task.get('status')

    def issue(self, sample, source, message):
        self.db.execute('INSERT INTO issues(sample_id,source,message) VALUES(?,?,?)', (sample, str(source), message))

    def source(self, path, role):
        # Reporting must never publish arbitrary files outside the chosen output tree.
        path.resolve().relative_to(self.root)
        relative = path.relative_to(self.root).as_posix()
        if relative in self.source_ids:
            return self.source_ids[relative]
        before = path.stat()
        self.source_stats[path] = (before.st_size, before.st_mtime_ns)
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        cur = self.db.execute('INSERT OR IGNORE INTO sources(path,bytes,sha256,role) VALUES(?,?,?,?)',
                             (relative, path.stat().st_size, digest.hexdigest(), role))
        identifier = self.db.execute('SELECT source_id FROM sources WHERE path=?', (relative,)).fetchone()[0]
        self.source_ids[relative] = identifier
        return identifier

    def one(self, sample, folder, pattern, role, required=False):
        matches = sorted(folder.glob(pattern))
        if len(matches) > 1:
            raise ValueError(f'Ambiguous {role} in {folder}: {[p.name for p in matches]}')
        if not matches:
            if required:
                self.issue(sample, folder.relative_to(self.root), f'Missing {role}; related fields may be unavailable.')
            return None
        self.source(matches[0], role)
        return matches[0]

    def orf(self, sample, identifier):
        if not identifier:
            raise ValueError(f'Empty ORF identifier in sample {sample}')
        self.db.execute('INSERT OR IGNORE INTO orfs(sample_id,orf_id) VALUES(?,?)', (sample, identifier))

    def sample(self, directory):
        s = directory.name if directory == self.root else directory.relative_to(self.root).as_posix()
        self.db.execute('INSERT INTO samples VALUES(?,?)', (s, directory.relative_to(self.root).as_posix()))
        self.db.execute('INSERT INTO entities(sample_id,entity_id,entity_type,pathway_status) VALUES(?,?,?,?)', (s, 'community', 'community', 'unavailable'))
        mapping = self.one(s, directory/'preprocessed', '*.mapping.txt', 'contig identifier mapping', True)
        if mapping:
            with mapping.open(newline='') as f:
                for r in csv.reader(f, delimiter='\t'):
                    if r:
                        if len(r) != 3:
                            raise ValueError(f'Expected three contig mapping columns: {mapping}')
                        self.db.execute('INSERT INTO contigs VALUES(?,?,?,?)', (s, r[0], r[1], number(r[2], True)))
        annotation_dir = directory/'results/annotation_table'
        primary = self.one(s, annotation_dir, '*.functional_and_taxonomic_table.txt', 'primary ORF annotations', True)
        if primary:
            for r in rows(primary):
                contig = r['Contig_Name']
                self.db.execute('INSERT OR IGNORE INTO contigs(sample_id,contig_id,length) VALUES(?,?,?)',
                                (s, contig, number(r['Contig_length'], True)))
                self.db.execute('INSERT INTO orfs VALUES(?,?,?,?,?,?,?,?,?,?,1)',
                    (s, r['ORF_ID'], contig, number(r['ORF_length'], True), number(r['start'], True),
                     number(r['end'], True), r['strand'], r['target'], r['product'], r['taxonomy']))
        # EC_RXN_map is the richer form of .1.txt; do not ingest both and duplicate hits.
        hits = self.one(s, annotation_dir, '*.EC_RXN_map.tsv', 'reference annotation and EC/reaction mapping')
        if hits is None:
            hits = self.one(s, annotation_dir, '*.1.txt', 'reference annotations', True)
        if hits:
            source_id = self.source(hits, 'reference annotations')
            last_orf = None
            for r in rows(hits):
                # MP's compact .1.txt leaves repeated ORF cells blank.
                last_orf = r['orf_id'] or last_orf
                self.orf(s, last_orf)
                cur = self.db.execute('INSERT INTO annotations(sample_id,orf_id,reference_db,target,product,score,ec,reaction,source_id) VALUES(?,?,?,?,?,?,?,?,?)',
                    (s, last_orf, r['ref dbname'], r['target'], r['product'], number(r['value']), r.get('EC',''), r.get('RXN',''), source_id))
                for kind, field in [('EC','EC'), ('reaction','RXN')]:
                    for term in r.get(field, '').split('|'):
                        if term and term != 'NONE':
                            self.db.execute('INSERT OR IGNORE INTO annotation_terms VALUES(?,?,?)', (cur.lastrowid, kind, term))
        groups = directory/'ptools/orf_map.txt'
        if groups.is_file():
            self.source(groups, 'collapsed ORF groups')
            with groups.open(newline='') as stream:
                for group in csv.reader(stream, delimiter='\t'):
                    group = [x for x in group if x]
                    if not group:
                        continue
                    for identifier in group:
                        self.orf(s, identifier)
                        self.db.execute('INSERT OR IGNORE INTO orf_groups VALUES(?,?,?)', (s, group[0], identifier))
        for mag in sorted((directory/'magsplitter/results').glob('*')):
            if not mag.is_dir():
                continue
            kind = 'unbinned' if 'non_binned' in mag.name else 'MAG'
            self.db.execute('INSERT OR IGNORE INTO entities(sample_id,entity_id,entity_type,pathway_status) VALUES(?,?,?,?)', (s, mag.name, kind, 'unavailable'))
            pf = mag/'0.pf'
            if pf.is_file():
                self.source(pf, 'MAG Pathway Tools input genes')
                with pf.open() as stream:
                    for line in stream:
                        if line.startswith('ID\t'):
                            identifier = line.rstrip('\r\n').split('\t', 1)[1]
                            self.orf(s, identifier)
                            self.db.execute('INSERT OR IGNORE INTO entity_orfs VALUES(?,?,?)', (s, mag.name, identifier))
        mag_map = directory/'magsplitter/contig_to_mag.tsv'
        if mag_map.is_file():
            source_id = self.source(mag_map, 'original contig-to-MAG membership')
            unmatched = 0
            with mag_map.open(newline='') as stream:
                for row in csv.reader(stream, delimiter='\t'):
                    if not row:
                        continue
                    if len(row) != 2 or not all(row):
                        raise ValueError(f'Expected original-contig and MAG columns without a header: {mag_map}')
                    original, original_mag_id = row
                    # MAGSplitter names output directories by replacing periods.
                    mag_id = original_mag_id.replace('.', '_')
                    self.db.execute('INSERT OR IGNORE INTO entities(sample_id,entity_id,entity_type,pathway_status) VALUES(?,?,?,?)', (s, mag_id, 'MAG', 'unavailable'))
                    contigs = self.db.execute('SELECT contig_id FROM contigs WHERE sample_id=? AND original_id=?', (s, original)).fetchall()
                    if not contigs:
                        unmatched += 1
                    for (contig,) in contigs:
                        self.db.execute('INSERT OR IGNORE INTO contig_mags VALUES(?,?,?,?,?)', (s, contig, mag_id, original_mag_id, source_id))
            if unmatched:
                self.issue(s, mag_map.relative_to(self.root), f'{unmatched} input contig assignments do not match retained contigs (for example after QC).')
        elif (directory/'magsplitter').is_dir():
            self.issue(s, 'magsplitter/contig_to_mag.tsv', 'Full contig-to-MAG membership is unavailable; MAG input genes alone are incomplete.')
        pgdb = directory/'results/pgdb'
        for path in sorted(pgdb.glob('community/*_pwy.tsv')) + sorted(pgdb.glob('MAGs/*/*_pwy.tsv')):
            entity = 'community' if path.parent.name == 'community' and path.parent.parent == pgdb else path.parent.name
            self.db.execute('INSERT OR IGNORE INTO entities(sample_id,entity_id,entity_type,pathway_status) VALUES(?,?,?,?)', (s, entity, 'community' if entity == 'community' else 'MAG', 'unavailable'))
            self.db.execute('UPDATE entities SET pathway_status=? WHERE sample_id=? AND entity_id=?', ('available', s, entity))
            status = self.pgdb_status.get((s, entity))
            if status and status not in ('SUCCESS', 'ALREADY_COMPUTED'):
                self.db.execute('UPDATE entities SET pathway_status=? WHERE sample_id=? AND entity_id=?', ('unavailable', s, entity))
                self.issue(s, path.relative_to(self.root), f'Pathway table excluded because the latest PGDB task is {status}.')
                continue
            source_id = self.source(path, 'pathway inference')
            for r in rows(path):
                self.db.execute('INSERT INTO pathways VALUES(?,?,?,?,?,?,?,?,?)',
                    (s, entity, r['PWY_NAME'], r['PWY_COMMON_NAME'], number(r['PWY_SCORE']),
                     number(r['NUM_REACTIONS'], True), number(r['NUM_COVERED_REACTIONS'], True), number(r['ORF_COUNT'], True), source_id))
                identifiers = {x.strip() for x in r['ORFS'].split(',') if x.strip() and x.strip() != 'nan'}
                for identifier in sorted(identifiers):
                    self.orf(s, identifier)
                    self.db.execute('INSERT INTO pathway_orfs VALUES(?,?,?,?)', (s, entity, r['PWY_NAME'], identifier))
                if number(r['ORF_COUNT'], True) != len(identifiers):
                    self.issue(s, path.relative_to(self.root), f"{entity}/{r['PWY_NAME']}: reported ORF count differs from unique listed ORFs.")
        for path in sorted((directory/'results/rpkm').glob('*.contig_counts.tsv')) + sorted((directory/'results/rpkm').glob('*.orf_counts.tsv')):
            feature = 'contig' if path.name.endswith('.contig_counts.tsv') else 'orf'
            source_id = self.source(path, 'read abundance')
            for r in rows(path):
                identifier = r['Contig'] if feature == 'contig' else r['Gene_ID']
                if feature == 'orf':
                    self.orf(s, identifier)
                    if r.get('seqname'):
                        self.db.execute('INSERT OR IGNORE INTO contigs(sample_id,contig_id) VALUES(?,?)', (s, r['seqname']))
                        self.db.execute('UPDATE orfs SET contig_id=COALESCE(contig_id,?) WHERE sample_id=? AND orf_id=?', (r['seqname'],s,identifier))
                else:
                    self.db.execute('INSERT OR IGNORE INTO contigs(sample_id,contig_id) VALUES(?,?)', (s,identifier))
                fields = {k:v for k,v in r.items() if k != 'Contig'} if feature == 'contig' else {k:r[k] for k in ('Count','RPKM','TPM')}
                for key, value in fields.items():
                    self.db.execute('INSERT INTO abundance VALUES(?,?,?,?,?,?)', (s, feature, identifier, key, number(value), source_id))
            self.issue(s, path.relative_to(self.root), 'Abundance is copied from source output; this report cannot establish that the original read mapping was correct.')
        missing = self.db.execute('SELECT COUNT(*) FROM orfs WHERE sample_id=? AND annotation_present=0', (s,)).fetchone()[0]
        if missing:
            self.issue(s, directory.relative_to(self.root), f'{missing} referenced ORF identifiers lack primary annotations; placeholders preserve these relationships.')
        self.db.commit()

    def inventory(self):
        for directory, dirs, files in os.walk(self.root, followlinks=False):
            dirs[:] = sorted(d for d in dirs if not d.startswith('.') and d not in ('reports','work','conda-cache') and not (Path(directory)/d).is_symlink())
            for name in sorted(files):
                path = Path(directory)/name
                if not path.is_file() or path.is_symlink() or name.startswith('.'):
                    continue
                stat = path.stat()
                self.db.execute('INSERT INTO files VALUES(?,?,?)',
                    (path.relative_to(self.root).as_posix(), stat.st_size, datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()))
        # Summaries record optional failures even when Nextflow completed normally.
        summaries = list(self.root.glob('logs/*/*/summary.json')) + list(self.root.glob('*/logs/*/*/summary.json'))
        for path in sorted(summaries, key=lambda p: p.stat().st_mtime_ns):
            source = path.relative_to(self.root).as_posix()
            data = json.loads(path.read_text())
            for task in data.get('tasks', []):
                sample_root = path.parent.parent.parent.parent
                sample_id = sample_root.name if sample_root == self.root else sample_root.relative_to(self.root).as_posix()
                if path.parent.parent.name == 'ptools' or task.get('entity'):
                    self.db.execute('UPDATE entities SET last_task_status=? WHERE sample_id=? AND entity_id=?',
                                    (task.get('status'), task.get('sample', sample_id), task.get('entity', task.get('label'))))
                self.db.execute('INSERT INTO execution VALUES(?,?,?,?,?,?,?,?)',
                    (path.parent.name, path.parent.parent.name, task.get('task',''), task.get('label',''),
                     task.get('status',''), task.get('elapsed_seconds'), task.get('error',''), source))
        self.db.commit()


def metadata(db):
    result = {}
    for name, (label, description) in VIEWS.items():
        result[name] = dict(label=label, description=description,
            columns=[dict(name=r[1], type=r[2]) for r in db.execute(f'PRAGMA table_info("{name}")')],
            rows=db.execute(f'SELECT COUNT(*) FROM "{name}"').fetchone()[0])
    tables = {}
    for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"):
        tables[name] = dict(columns=[dict(name=r[1], type=r[2], primary_key_order=r[5]) for r in db.execute(f'PRAGMA table_info("{name}")')],
                           foreign_keys=[dict(id=r[0], sequence=r[1], table=r[2], column=r[3], references=r[4]) for r in db.execute(f'PRAGMA foreign_key_list("{name}")')])
    return dict(schema_version=SCHEMA_VERSION, generated_utc=datetime.now(timezone.utc).isoformat(), views=result, tables=tables)


def atomic_text(path, value):
    temporary = path.with_name('.'+path.name+'.tmp')
    temporary.write_text(value, encoding='utf-8')
    temporary.replace(path)


def write_report(db, root, reports, info):
    esc = html.escape
    cards = ''.join(f'<li><a href="EDA_portal.html#table={quote(name)}">{esc(view["label"])}</a>: {view["rows"]:,} rows</li>' for name, view in info['views'].items())
    samples = ''.join(f'<tr><td>{esc(s)}</td><td>{esc(p)}</td></tr>' for s,p in db.execute('SELECT * FROM samples ORDER BY sample_id'))
    notes = ''.join(f'<li>{esc(str(s))}: {esc(m)}</li>' for s,m in db.execute('SELECT sample_id,message FROM issues LIMIT 100'))
    # Link files rather than relying on web-server directory listings.
    links = ''.join(f'<li><a href="../{quote(p, safe="/")}">{esc(p)}</a></li>' for (p,) in db.execute("SELECT path FROM files WHERE path LIKE '%/trace.tsv' OR path LIKE '%/report.html' OR path LIKE '%/timeline.html' OR path LIKE '%/summary.json' ORDER BY path"))
    document = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>MetaPathways run report</title>
<style>body{{font:17px system-ui,sans-serif;max-width:1100px;margin:3rem auto;padding:0 1rem;color:#203038}}a{{color:#006a80}}table{{border-collapse:collapse}}th,td{{padding:.6rem;border-bottom:1px solid #ddd;text-align:left}}li{{margin:.5rem 0}}code{{background:#eef4f5;padding:.15rem}}</style></head>
<body><h1>MetaPathways run report</h1><p>Generated {esc(info['generated_utc'])}. An inventory of existing outputs; it does not rerun or interpret analyses.</p>
<p><a href="EDA_portal.html">Open EDA portal</a> · <a href="results.sqlite">Relational results database</a> · <a href="schema.json">Schema and table definitions</a> · <a href="output_inventory.tsv">All output files</a></p>
<p>For searching and CSV export: <code>metapathways report -o OUTPUT --serve --no-rebuild</code>. The HTML report and file links also work without the server.</p>
<h2>Samples</h2><table><tr><th>Sample</th><th>Output directory</th></tr>{samples}</table>
<h2>Available results</h2><ul>{cards}</ul><p>Counts describe table rows, not necessarily distinct genes or pathways. Sample and entity identifiers scope every biological relationship.</p>
<h2>Availability and import notes</h2><ul>{notes or '<li>No import notes.</li>'}</ul><p>Up to 100 notes shown; the portal includes all notes. Missing pathways are not evidence of biological absence. Expected MAG failures remain in execution history.</p>
<h2>Nextflow run details</h2><ul>{links or '<li>No Nextflow records found in these outputs.</li>'}</ul></body></html>'''
    atomic_text(reports/'MP_run_report.html', document)
    with (reports/'output_inventory.tsv').open('w', newline='') as f:
        writer=csv.writer(f, delimiter='\t'); writer.writerow(['path','bytes','modified_utc'])
        writer.writerows(db.execute('SELECT * FROM files ORDER BY path'))


def build_report(output):
    root = Path(output).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f'Output directory does not exist: {root}')
    samples = sample_dirs(root)
    if not samples:
        raise ValueError(f'No sample outputs found directly in {root} or its immediate child directories')
    reports = root/'reports'
    reports.mkdir(exist_ok=True)
    with (reports/'.build.lock').open('w') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Another report is being built in {reports}')
        with tempfile.TemporaryDirectory(prefix='.build-', dir=reports) as temporary:
            database = Path(temporary)/'results.sqlite'
            with closing(sqlite3.connect(database)) as db:
                # Build only in a disposable staging file; publish after validation.
                db.execute('PRAGMA journal_mode=MEMORY')
                db.execute('PRAGMA synchronous=OFF')
                db.execute('PRAGMA foreign_keys=ON')
                db.executescript('BEGIN;\n' + SCHEMA)
                db.execute('PRAGMA cache_size=-32768')
                importer = Importer(db, root)
                for directory in samples:
                    print(f'Indexing report tables: {directory.name}', flush=True)
                    importer.sample(directory)
                importer.inventory()
                for source, expected in importer.source_stats.items():
                    current = source.stat()
                    if (current.st_size, current.st_mtime_ns) != expected:
                        raise ValueError(f"Source changed during report generation: {source}; rebuild when outputs are stable")
                if db.execute('PRAGMA foreign_key_check').fetchone():
                    raise ValueError('Report database contains invalid relationships')
                info = metadata(db)
                info['output_root'] = str(root)
                info['sample_paths'] = [str(p) for p in samples]
                db.execute('ANALYZE'); db.commit()
                write_report(db, root, reports, info)
            with database.open('rb') as stream:
                os.fsync(stream.fileno())
            database.replace(reports/'results.sqlite')
            atomic_text(reports/'schema.json', json.dumps(info, indent=2)+'\n')
            assets = Path(__file__).parent/'report_assets'
            for name in ('EDA_portal.html','portal.js','portal.css'):
                shutil.copyfile(assets/name, reports/name)
    print(f'Report: {reports / "MP_run_report.html"}', flush=True)
    return reports
