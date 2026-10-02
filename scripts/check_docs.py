#!/usr/bin/env python3
"""Check repository documentation links and the generated CLI reference."""
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT=Path(__file__).resolve().parents[1]
files=[ROOT/'README.md', *sorted((ROOT/'docs').glob('*.md')), ROOT/'docker/README.quay.md']
errors=[]
for source in files:
    content=re.sub(r'```.*?```','',source.read_text(),flags=re.S)
    for link in re.findall(r'\]\(([^\s)]+)\)',content):
        url=urlsplit(link)
        if url.scheme or url.netloc or not url.path:
            continue
        target=(source.parent/unquote(url.path)).resolve()
        if not target.exists():
            errors.append(f'{source.relative_to(ROOT)}: missing link {link}')
if errors:
    raise SystemExit('\n'.join(errors))
subprocess.run([sys.executable,str(ROOT/'scripts/generate_cli_docs.py'),'--check'],check=True)
print(f'Local links checked in {len(files)} Markdown documents.')
