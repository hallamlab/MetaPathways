"""Loopback-only, read-only SQL queries and streamed CSV exports for MP reports."""
import argparse
from contextlib import closing
import csv
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
import json
from pathlib import Path
import secrets
import sqlite3
import threading
import time
from urllib.parse import parse_qs, unquote, urlsplit, quote
import webbrowser

from metapathways.reporting import VIEWS, build_report


def identifier(value):
    return '"' + value.replace('"', '""') + '"'


def connect(database):
    db = sqlite3.connect(database.as_uri() + '?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA query_only=ON')
    deadline = time.monotonic() + 120
    db.set_progress_handler(lambda: int(time.monotonic() > deadline), 10000)
    return db


def query(db, spec, export=False):
    table = spec.get('table', 'orf_explorer')
    if table not in VIEWS:
        raise ValueError('Unknown table')
    details = list(db.execute(f'PRAGMA table_info({identifier(table)})'))
    columns = [r[1] for r in details]
    numeric = {r[1] for r in details if r[2] in ('REAL','INTEGER')} | {'linked_orf_count'}
    selected = spec.get('columns') or columns
    if not isinstance(selected, list) or len(selected) > len(columns) or any(c not in columns for c in selected):
        raise ValueError('Unknown selected column')
    where, values = [], []
    search = spec.get('search', '')
    if not isinstance(search, str) or len(search) > 1000:
        raise ValueError('Search is too long')
    if search:
        where.append('(' + ' OR '.join(f'instr(lower(CAST(t.{identifier(c)} AS TEXT)),lower(?))>0' for c in columns) + ')')
        values.extend([search] * len(columns))
    filters = spec.get('filters', [])
    if not isinstance(filters, list) or len(filters) > 30:
        raise ValueError('At most 30 filters are supported')
    for f in filters:
        column, op, value = f.get('column'), f.get('op'), f.get('value', '')
        if column not in columns:
            raise ValueError('Unknown filter column')
        expr = 't.' + identifier(column)
        if op == 'missing':
            where.append(f'({expr} IS NULL OR {expr}=\'\')'); continue
        if op == 'present':
            where.append(f'({expr} IS NOT NULL AND {expr}!=\'\')'); continue
        if not isinstance(value, (str, int, float)) or len(str(value)) > 1000:
            raise ValueError('Invalid filter value')
        if op == 'contains':
            where.append(f'instr(lower(CAST({expr} AS TEXT)),lower(?))>0')
        elif op in ('eq','ne'):
            where.append(f'{expr} {"=" if op == "eq" else "!="} ?')
        elif op in ('lt','le','gt','ge'):
            if column not in numeric:
                raise ValueError('Numeric comparisons require a numeric column')
            try:
                value = float(value)
            except (ValueError, TypeError):
                raise ValueError('Numeric comparison requires a number')
            where.append(f'CAST({expr} AS REAL) ' + {'lt':'<','le':'<=','gt':'>','ge':'>='}[op] + ' ?')
        else:
            raise ValueError('Unknown filter operator')
        values.append(value)
    # Related filters use EXISTS: keep one ORF/annotation row rather than multiply it.
    related = spec.get('related') or {}
    if not isinstance(related, dict) or set(related) - {'pathway_id','entity_id','ec','reaction','reference_db'}:
        raise ValueError('Unknown related filter')
    if related:
        if not {'sample_id','orf_id'} <= set(columns):
            raise ValueError('Related filters require an ORF-based table')
        if any(not isinstance(v, str) or len(v) > 1000 for v in related.values()):
            raise ValueError('Invalid related filter')
        if related.get('pathway_id'):
            expression = 'EXISTS(SELECT 1 FROM pathway_orfs p WHERE p.sample_id=t.sample_id AND p.orf_id=t.orf_id AND p.pathway_id=?'
            values.append(related['pathway_id'])
            if related.get('entity_id'):
                expression += ' AND p.entity_id=?'; values.append(related['entity_id'])
            where.append(expression + ')')
        elif related.get('entity_id'):
            where.append('EXISTS(SELECT 1 FROM entity_orfs e WHERE e.sample_id=t.sample_id AND e.orf_id=t.orf_id AND e.entity_id=?)')
            values.append(related['entity_id'])
        if any(related.get(k) for k in ('ec','reaction','reference_db')):
            expression = 'EXISTS(SELECT 1 FROM annotations a WHERE a.sample_id=t.sample_id AND a.orf_id=t.orf_id'
            if related.get('reference_db'):
                expression += ' AND a.reference_db=?'; values.append(related['reference_db'])
            for kind in ('ec','reaction'):
                if related.get(kind):
                    expression += ' AND EXISTS(SELECT 1 FROM annotation_terms x WHERE x.annotation_id=a.annotation_id AND x.term_type=? AND x.term=?)'
                    values.extend(['EC' if kind == 'ec' else 'reaction', related[kind]])
            where.append(expression + ')')
    base = ' FROM ' + identifier(table) + ' t' + (' WHERE ' + ' AND '.join(where) if where else '')
    sort = spec.get('sort') or columns[0]
    if sort not in columns:
        raise ValueError('Unknown sort column')
    order = ' ORDER BY t.' + identifier(sort) + (' DESC' if spec.get('descending') else ' ASC')
    for column in ('sample_id','orf_id','entity_id','pathway_id','annotation_id','source_id'):
        if column in columns and column != sort:
            order += ',t.' + identifier(column)
    sql = 'SELECT ' + ','.join('t.'+identifier(c) for c in selected) + base + order
    if export:
        return selected, db.execute(sql, values)
    limit, offset = int(spec.get('limit',100)), int(spec.get('offset',0))
    if not 1 <= limit <= 500 or offset < 0:
        raise ValueError('Invalid page limits')
    total = db.execute('SELECT COUNT(*)' + base, values).fetchone()[0]
    result = db.execute(sql + ' LIMIT ? OFFSET ?', values+[limit,offset])
    return dict(columns=selected, rows=[dict(r) for r in result], total=total, offset=offset, limit=limit)


def csv_value(value):
    # Keep numeric negatives numeric, and prevent text from becoming spreadsheet formulas.
    if isinstance(value, str) and value.lstrip().startswith(('=', '+', '-', '@')):
        return "'" + value
    return value


class ReportServer(ThreadingHTTPServer):
    daemon_threads = True
    def __init__(self, root, port=0):
        self.root = Path(root).resolve()
        self.token = secrets.token_urlsafe(24)
        self.slots = threading.BoundedSemaphore(4)
        super().__init__(('127.0.0.1', port), Handler)

    @property
    def url(self):
        return f'http://127.0.0.1:{self.server_port}/{self.token}/reports/EDA_portal.html'


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        # Queries can include scientific identifiers; keep the terminal concise.
        if len(args) > 1 and str(args[1]).startswith(('4','5')):
            print('Portal request failed:', args[1], flush=True)

    def send(self, value, kind='application/json', status=200):
        raw = json.dumps(value).encode() if kind == 'application/json' else value
        self.send_response(status)
        self.send_header('Content-Type', kind)
        self.send_header('Content-Length', str(len(raw)))
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        self.end_headers(); self.wfile.write(raw)

    def do_GET(self):
        try:
            self.handle_get()
        except (BrokenPipeError, ConnectionResetError):
            pass
        except (ValueError, KeyError, TypeError, sqlite3.Error, OSError) as exc:
            if getattr(self, '_response_started', False):
                self.close_connection = True
            else:
                self.send(dict(error=str(exc)), status=400)

    def send_response(self, code, message=None):
        self._response_started = True
        super().send_response(code, message)

    def handle_get(self):
        expected = f'127.0.0.1:{self.server.server_port}'
        if self.headers.get('Host') != expected or self.headers.get('Origin') not in (None, f'http://{expected}'):
            self.send(dict(error='Local same-origin access only'), status=403); return
        parsed = urlsplit(self.path)
        prefix = '/' + self.server.token + '/'
        if not parsed.path.startswith(prefix):
            self.send(dict(error='Not found'), status=404); return
        path = unquote(parsed.path[len(prefix):])
        root = self.server.root
        if path == 'api/meta':
            self.send(json.loads((root/'reports/schema.json').read_text())); return
        if path in ('api/query','api/export'):
            raw = parse_qs(parsed.query).get('spec',['{}'])[0]
            if len(raw) > 20000:
                raise ValueError('Query is too large')
            spec = json.loads(raw)
            if not isinstance(spec, dict):
                raise ValueError('Query must be an object')
            if not self.server.slots.acquire(blocking=False):
                self.send(dict(error='Four queries are already active; try again shortly'), status=429); return
            try:
                with closing(connect(root/'reports/results.sqlite')) as db:
                    if path == 'api/query':
                        self.send(query(db, spec)); return
                    columns, cursor = query(db, spec, export=True)
                    self.send_response(200)
                    self.send_header('Content-Type','text/csv; charset=utf-8')
                    self.send_header('Content-Disposition','attachment; filename="metapathways-subset.csv"')
                    self.send_header('Cache-Control','no-store')
                    self.end_headers()
                    buffer=io.StringIO(); writer=csv.writer(buffer)
                    writer.writerow(columns)
                    for row in cursor:
                        writer.writerow([csv_value(v) for v in row])
                        if buffer.tell() > 65536:
                            self.wfile.write(buffer.getvalue().encode()); buffer.seek(0); buffer.truncate()
                    self.wfile.write(buffer.getvalue().encode())
            finally:
                self.server.slots.release()
            return
        target=(root/path).resolve()
        try:
            relative=target.relative_to(root)
        except ValueError:
            self.send(dict(error='Not found'), status=404); return
        if any(p.startswith('.') for p in relative.parts) or not target.is_file():
            self.send(dict(error='Not found'), status=404); return
        if relative.parts[0] != 'reports':
            with closing(connect(root/'reports/results.sqlite')) as db:
                if not db.execute('SELECT 1 FROM files WHERE path=?',(relative.as_posix(),)).fetchone():
                    self.send(dict(error='Not inventoried'), status=404); return
        import mimetypes
        self.send_response(200)
        self.send_header('Content-Type',mimetypes.guess_type(target.name)[0] or 'application/octet-stream')
        self.send_header('Content-Length',str(target.stat().st_size))
        self.send_header('X-Content-Type-Options','nosniff')
        self.end_headers()
        with target.open('rb') as stream:
            for block in iter(lambda: stream.read(65536), b''):
                self.wfile.write(block)


def main(argv=None):
    parser=argparse.ArgumentParser(description='Build a navigable report from existing MP outputs; no analyses are rerun.')
    parser.add_argument('-o','--output_dir',required=True,help='sample output or parent containing sample outputs')
    parser.add_argument('--serve',action='store_true',help='open the local searchable portal and serve until Ctrl-C')
    parser.add_argument('--no-rebuild',action='store_true',help='use an existing report snapshot')
    parser.add_argument('--port',type=int,default=0,help='local port [automatically selected]')
    parser.add_argument('--no-browser',action='store_true',help='print the local URL without opening a browser')
    args=parser.parse_args(argv)
    root=Path(args.output_dir).expanduser().resolve()
    if not args.no_rebuild:
        build_report(root)
    if not (root/'reports/results.sqlite').is_file():
        parser.error('No report database exists; omit --no-rebuild')
    if args.serve:
        with ReportServer(root,args.port) as server:
            print(f'EDA portal: {server.url}\nLocal access only. Ctrl-C stops the portal.',flush=True)
            if not args.no_browser:
                webbrowser.open(server.url)
            try:
                server.serve_forever()
            except KeyboardInterrupt:
                print('Portal stopped.',flush=True)
