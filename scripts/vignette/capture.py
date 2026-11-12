"""Capture the actual example explorer and clearly labelled command references.

Requires Playwright and a built example. No browser address/token is captured.
"""
import argparse,csv,json,sys,threading
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from metapathways.report_server import ReportServer
from playwright.sync_api import sync_playwright
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('example',type=Path)
parser.add_argument('--browser',required=True)
args=parser.parse_args()
assets=Path(__file__).resolve().parents[2]/'docs/assets/vignette'
server=ReportServer(args.example/'all')
threading.Thread(target=server.serve_forever,daemon=True).start()
exports=args.example/'exports';exports.mkdir(exist_ok=True)
try:
 with sync_playwright() as p:
  browser=p.chromium.launch(executable_path=args.browser,args=['--no-sandbox'])
  page=browser.new_page(viewport={'width':1560,'height':1050},device_scale_factor=1)
  errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
  page.goto(server.url);page.wait_for_selector('#table option',state='attached')
  def shot(name):page.screenshot(path=str(assets/(name+'.png')),full_page=True)
  def query(table,filters=(),columns=None,search=''):
   page.evaluate('''async ({table,filters,columns,search})=>{configure(table,filters);document.getElementById('search').value=search;if(columns){document.querySelectorAll('#columns input').forEach(c=>c.checked=columns.includes(c.value));}await apply();}''',dict(table=table,filters=list(filters),columns=columns,search=search))
   page.wait_for_function("document.querySelector('#status').textContent.includes('matching rows')")
  def eq(column,value):return dict(column=column,op='eq',value=value)
  def export(name,count):
   with page.expect_download() as d:page.locator('#export').click()
   d.value.save_as(str(exports/name))
   rows=list(csv.DictReader((exports/name).open()));assert len(rows)==count,(name,len(rows));return rows
  shot('01-samples')
  query('annotation_explorer',columns=['sample_id','orf_id','contig_id','reference_db','product','ec','reaction'],search='pantothenate synthetase')
  assert page.locator('#status').inner_text()=='2 matching rows'
  shot('02-function');export('panc-annotations.csv',2)
  sample='Gastrointestinal_5';contig='Gastrointestinal_5-C2';mag='CAMI_OTU_97_68_0'
  query('contig_mags',[eq('sample_id',sample),eq('contig_id',contig)],['sample_id','contig_id','entity_id','original_mag_id']);shot('03-membership')
  query('rna_taxonomy_explorer',[eq('sample_id',sample),eq('contig_id',contig)],['sample_id','contig_id','subunit','similarity','taxonomy']);shot('04-rna')
  query('pathway_explorer',[eq('sample_id',sample),eq('entity_id',mag)],['sample_id','entity_id','pathway_id','name','reactions','covered_reactions','linked_orf_count'],search='PANTO-PWY');shot('05-pathway')
  query('pathway_gene_explorer',[eq('sample_id',sample),eq('entity_id',mag),eq('pathway_id','PANTO-PWY')],['sample_id','entity_id','pathway_id','orf_id','contig_id','product']);shot('06-genes');export('panto-pathway-genes.csv',2)
  query('orf_explorer',[eq('sample_id',sample),eq('contig_id',contig),dict(column='start',op='ge',value='10576'),dict(column='start',op='le',value='12296')],['sample_id','orf_id','start','end','product']);shot('07-neighborhood')
  with page.expect_popup() as popup:page.get_by_role('button',name='FASTA',exact=True).first.click()
  fasta=popup.value;fasta.wait_for_selector('#record option:nth-child(2)',state='attached');fasta.select_option('#record',index=1);fasta.wait_for_function("document.querySelector('#fasta').textContent.startsWith('>')");fasta.screenshot(path=str(assets/'08-fasta.png'),full_page=True);fasta.close()
  page.goto(server.url.replace('EDA_portal.html','downloads.html'));page.wait_for_selector('#selections input')
  page.locator(f'input[data-sample="{sample}"][data-mag="{mag}"]').check()
  page.locator('input[name="mode"][value="complete"]').check();page.locator('#preview').click();page.wait_for_function("!document.getElementById('download').disabled")
  shot('09-download')
  with page.expect_download() as d:page.locator('#download').click()
  d.value.save_as(str(exports/'roseburia-fragment.zip'))
  import zipfile
  with zipfile.ZipFile(exports/'roseburia-fragment.zip') as z:
   assert any(n.endswith('.gbk') for n in z.namelist()),z.namelist()
   assert any(n.endswith('.fasta') for n in z.namelist())
  # Command references contain no fabricated output. Historical/verified output
  # is separately labelled; copyable equivalents live in the vignette.
  commands=json.loads((assets/'terminal-steps.json').read_text())
  for step in commands:
   page.set_content('<html><body style="margin:0;background:#16202b;color:#f3f5f7;font:18px monospace;padding:28px"><div id="label" style="color:#a8d5ef;margin-bottom:24px"></div><pre id="terminal" style="white-space:pre-wrap;line-height:1.5"></pre></body></html>')
   page.locator('#label').evaluate('(e,t)=>e.textContent=t',step['label'])
   page.locator('#terminal').evaluate('(e,t)=>e.textContent=t',step['text'])
   page.set_viewport_size({'width':1200,'height':max(220,110+28*len(step['text'].splitlines()))})
   page.screenshot(path=str(assets/(step['id']+'.png')),full_page=True)
  assert not errors,errors
  browser.close()
  print('Verified PanC annotations, pathway exports, relocated explorer, and MAG ZIP with GenBank/FASTA. Screenshots saved.')
finally:server.shutdown();server.server_close()
