"""Opt-in Chromium regression for explorer scope, keyboard filters and reset layout."""
import argparse
import sys,threading,csv
from pathlib import Path
root=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(root),str(root/'tests')]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--browser-executable', default=None)
args=parser.parse_args()
from test_report_downloads import DownloadTests
from metapathways.report_server import ReportServer
from playwright.sync_api import sync_playwright
fixture=DownloadTests();fixture.setUp()
server=ReportServer(fixture.root);threading.Thread(target=server.serve_forever,daemon=True).start()
try:
    with sync_playwright() as p:
        browser=p.chromium.launch(executable_path=args.browser_executable,args=['--no-sandbox']);page=browser.new_page(viewport={'width':1400,'height':1000})
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(server.url)
        page.wait_for_selector('#table option',state='attached')
        page.select_option('#scopeSample','alpha');page.select_option('#scopeMag','MAG1')
        page.wait_for_function("document.querySelector('#status').textContent==='1 matching rows'")
        assert page.locator('#scopeSample').input_value()=='alpha'
        page.select_option('#table','entities')
        page.wait_for_function("document.querySelector('#status').textContent==='1 matching rows'")
        page.select_option('#table','files')
        page.wait_for_function("document.querySelector('#status').textContent.includes('matching rows')")
        page.reload()
        page.wait_for_function("document.querySelector('#scopeMag').value==='MAG1' && document.querySelector('#status').textContent.includes('matching rows')")
        assert page.locator('#scopeSample').input_value()=='alpha'
        page.select_option('#table','mag_orf_explorer')
        page.wait_for_function("document.querySelector('#status').textContent==='2 matching rows'")
        with page.expect_download() as download:
            page.get_by_role('button',name='Export all matching rows as CSV',exact=True).click()
        download.value.save_as(str(fixture.root/'scope.csv'))
        with open(str(fixture.root/'scope.csv')) as f:
            rows=list(csv.DictReader(f))
        assert len(rows)==2 and {r['sample_id'] for r in rows}=={'alpha'} and {r['entity_id'] for r in rows}=={'MAG1'}
        page.select_option('#table','orf_explorer')
        page.get_by_role('button',name='Add column filter',exact=True).click()
        page.get_by_label('Filter column',exact=True).select_option('product')
        page.get_by_label('Filter value',exact=True).fill('kin')
        page.get_by_label('Filter value',exact=True).press('Enter')
        page.wait_for_function("document.querySelector('#status').textContent==='1 matching rows'")
        page.locator('#search').fill('no match');page.locator('#search').press('Enter')
        page.wait_for_function("document.querySelector('#status').textContent==='0 matching rows'")
        page.locator('#search').fill('');page.locator('#search').press('Enter')
        page.wait_for_function("document.querySelector('#status').textContent==='1 matching rows'")
        page.locator('#relatedPanel summary').click()
        page.locator('#related_entity_id').fill('missing');page.locator('#related_entity_id').press('Enter')
        page.wait_for_function("document.querySelector('#status').textContent==='0 matching rows'")
        page.locator('#related_entity_id').fill('MAG1');page.locator('#related_entity_id').press('Enter')
        page.wait_for_function("document.querySelector('#status').textContent==='1 matching rows'")
        sizes=page.evaluate("['resetScope','reset'].map(id=>{const r=document.getElementById(id).getBoundingClientRect();return [r.width,r.height]})")
        assert all(abs(a-b)<.1 for a,b in zip(*sizes)),sizes
        assert not errors,errors
        browser.close()
        print('Passed: repeated MAG IDs isolated by sample in browsing/CSV, scope survives metadata/reload, Enter applies all three search types, both reset buttons match dimensions.')
finally:
    server.shutdown();server.server_close();fixture.doCleanups()
