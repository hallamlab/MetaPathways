#!/usr/bin/env python3
"""Rebuild shared-scale SVG previews from Mermaid and publication SVG sources.

Run from the repository root after installing docs/diagram-requirements.txt
and running `python -m playwright install chromium`. Rendering downloads the
pinned Mermaid bundle; documentation builds use the committed SVGs offline.
"""
from pathlib import Path
import json
import argparse
import xml.etree.ElementTree as ET
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
NS = '{http://www.w3.org/2000/svg}'
ET.register_namespace('', NS[1:-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mermaid-js", type=Path, help="Local Mermaid 11.12.1 bundle")
    parser.add_argument("--supporting-only", action="store_true", help="Rebuild Mermaid diagrams without changing the primary figures")
    args = parser.parse_args()
    config = json.loads((ROOT / 'docs/diagrams/figures.json').read_text())
    from render_detailed_workflow import render, build
    if not args.supporting_only:
        build(ROOT)
        from render_nodal_brief import build as build_nodal
        import cairosvg
        nodal = build_nodal()
        cairosvg.svg2pdf(url=str(nodal),write_to=str(nodal.with_suffix(".pdf")))
        detailed_svg, detailed_width, detailed_height = render(ROOT)
    rendered = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto('about:blank')
        if args.mermaid_js:
            page.add_script_tag(path=str(args.mermaid_js.resolve()))
        else:
            page.add_script_tag(url='https://cdn.jsdelivr.net/npm/mermaid@11.12.1/dist/mermaid.min.js')
        page.evaluate('mermaid.initialize({startOnLoad:false,securityLevel:"strict",htmlLabels:false,flowchart:{htmlLabels:false}})')
        for item in config['figures']:
            code = (ROOT / 'docs/diagrams' / (item['name'] + '.mmd')).read_text()
            svg = page.evaluate('''async ({code,id}) => {
                const result = await mermaid.render(id, code);
                const div = document.createElement('div');
                div.innerHTML = result.svg;
                return new XMLSerializer().serializeToString(div.firstElementChild);
            }''', {'code': code, 'id': ROOT.name.lower() + '-' + item['name']})
            element = ET.fromstring(svg)
            if list(element.iter(NS + 'foreignObject')):
                raise RuntimeError(f"{item['name']}: HTML labels are not portable; require native SVG text")
            for group in element.iter(NS + 'g'):
                if 'node' in group.get('class', '').split() and not any(
                    ''.join(t.itertext()).strip() for t in group.iter(NS + 'text')
                ):
                    raise RuntimeError(f"{item['name']}: node lacks a native text label")
            rendered.append((item['name'], element, config['mermaid_scale']))
        if not args.supporting_only:
            page.set_content('<html><head><style>html,body{margin:0;padding:0}svg{display:block}</style></head><body>' + detailed_svg.read_text() + '</body></html>')
            page.pdf(path=str(detailed_svg.with_suffix('.pdf')), width=f'{detailed_width}px', height=f'{detailed_height}px', print_background=True, margin=dict(top='0', right='0', bottom='0', left='0'))
            brief = ROOT / 'docs/assets/workflow-main.svg'
            page.set_content('<html><head><style>html,body{margin:0;padding:0}svg{display:block}</style></head><body>' + brief.read_text() + '</body></html>')
            page.pdf(path=str(brief.with_suffix('.pdf')), width=ET.parse(brief).getroot().get('width') + 'px', height=ET.parse(brief).getroot().get('height') + 'px', print_background=True, margin=dict(top='0', right='0', bottom='0', left='0'))
        browser.close()
    for name in ([] if args.supporting_only else ['workflow', 'workflow-brief', 'workflow-detailed']):
        source_name = 'workflow-main' if name == 'workflow' else name
        source = ROOT / 'docs/assets' / (source_name + '.svg')
        if source.exists():
            rendered.append((name, ET.parse(source).getroot(), 1))
    width = config['canvas_width']
    destination = ROOT / 'docs/assets/diagrams'
    destination.mkdir(exist_ok=True)
    # Check every diagram before replacing any previews. Keep the common width
    # synchronized between ASPIRE and MetaPathways if a future figure needs more.
    for name, element, scale in rendered:
        if name in ('workflow', 'workflow-detailed'):
            continue  # Primary figures preserve their own complete canvas.
        natural_width = float(element.get('viewBox').split()[2]) * scale
        if natural_width > width - 40:
            raise SystemExit(f'{name} needs {natural_width + 40:g} canvas units; update the shared canvas width.')
    for name, element, scale in rendered:
        if name in ('workflow', 'workflow-detailed'):
            # Primary figures match ASPIRE v1: preserve the complete source canvas.
            source_name = 'workflow-main' if name == 'workflow' else name
            (destination / (name + '.svg')).write_text((ROOT / 'docs/assets' / (source_name + '.svg')).read_text())
            continue
        _, _, w, h = map(float, element.get('viewBox').split())
        w *= scale
        h *= scale
        canvas = ET.Element(NS + 'svg', {
            'viewBox': f'0 0 {width} {h + 40:g}', 'width': str(width),
            'height': f'{h + 40:g}', 'role': 'img',
        })
        ET.SubElement(canvas, NS + 'title').text = name.replace('-', ' ').title()
        ET.SubElement(canvas, NS + 'rect', {
            'width': str(width), 'height': f'{h + 40:g}', 'fill': '#ffffff',
        })
        element.set('x', f'{(width - w) / 2:g}')
        element.set('y', '20')
        element.set('width', f'{w:g}')
        element.set('height', f'{h:g}')
        element.attrib.pop('style', None)
        canvas.append(element)
        (destination / (name + '.svg')).write_text(ET.tostring(canvas, encoding='unicode') + '\n')
        print(f'Rendered {name} on {width}-unit canvas')


if __name__ == '__main__':
    main()
