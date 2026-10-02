#!/usr/bin/env python3
"""Validate local links and fragments in the built documentation."""
import argparse
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        if tag == 'a' and 'href' in attrs:
            self.links.append(attrs['href'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--readme', type=Path, help='also check hosted documentation links in the README')
    args = parser.parse_args()
    root = args.directory.resolve()
    pages = {p.resolve(): Page(p.read_text()) for p in root.rglob('*.html')}
    if not pages:
        raise SystemExit(f'No HTML pages found in {root}')
    errors = []
    for source, page in pages.items():
        for link in page.links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            target = (source.parent / unquote(url.path)).resolve() if url.path else source
            if target.is_dir():
                target /= 'index.html'
            if not target.exists():
                errors.append(f'{source.name}: missing {link}')
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                errors.append(f'{source.name}: missing anchor {link}')
    if args.readme:
        for link in re.findall(r'\]\((https://metapathways\.readthedocs\.io/[^\s)]*)\)', args.readme.read_text()):
            url = urlsplit(link)
            relative = url.path.removeprefix('/en/latest/').lstrip('/') or 'index.html'
            target = root / relative
            if target not in pages:
                errors.append(f'README: missing documentation page {link}')
            elif url.fragment and unquote(url.fragment) not in pages[target].ids:
                errors.append(f'README: missing documentation anchor {link}')
    if errors:
        raise SystemExit('\n'.join(sorted(set(errors))))
    print(f'Local links and anchors checked in {len(pages)} HTML pages.')


if __name__ == '__main__':
    main()
