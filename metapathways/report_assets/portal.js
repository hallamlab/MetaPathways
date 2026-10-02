'use strict';
const $ = id => document.getElementById(id);
let metadata, active, offset = 0, total = 0, sort = '', descending = false, serial = 0;
const relatedKeys = ['pathway_id','entity_id','ec','reaction','reference_db'];
function element(tag, text) { const e=document.createElement(tag); if(text!==undefined)e.textContent=text; return e; }
function columns() { return metadata.views[$('table').value].columns.map(c=>c.name); }
function button(text, action) { const b=element('button',text); b.type='button'; b.onclick=action; return b; }
function api(path) { return new URL('../api/'+path, location.href).href; }
function queryURL(kind, spec) { return api(kind)+'?spec='+encodeURIComponent(JSON.stringify(spec)); }
function option(value, label=value) { const o=element('option',label); o.value=value; return o; }
function filterRow(column=columns()[0],op='contains',value='') {
  const row=element('div'); row.className='filter';
  const c=element('select'); c.setAttribute('aria-label','Filter column'); columns().forEach(x=>c.append(option(x))); c.value=column;
  const operator=element('select'); operator.setAttribute('aria-label','Comparison');
  Object.entries({contains:'contains',eq:'equals',ne:'does not equal',ge:'≥',le:'≤',gt:'>',lt:'<',missing:'is missing',present:'is present'}).forEach(([k,v])=>operator.append(option(k,v))); operator.value=op;
  const input=element('input'); input.value=value; input.setAttribute('aria-label','Filter value');
  row.append(c,operator,input,button('Remove',()=>row.remove())); $('filters').append(row);
}
function configure(table, filters=[]) {
  if(!metadata.views[table])table='samples';
  $('table').value=table; $('description').textContent=metadata.views[table].description;
  $('filters').replaceChildren(); $('columns').replaceChildren(); $('search').value='';
  relatedKeys.forEach(k=>$('related_'+k).value='');
  $('relatedPanel').hidden=!(columns().includes('sample_id')&&columns().includes('orf_id'));
  columns().forEach(c=>{const label=element('label');const box=element('input');box.type='checkbox';box.value=c;box.checked=true;label.append(box,document.createTextNode(c));$('columns').append(label);});
  filters.forEach(f=>filterRow(f.column,f.op,f.value)); sort='';descending=false;offset=0;
}
function specFromControls() {
  const filters=Array.from($('filters').children).map(row=>({column:row.children[0].value,op:row.children[1].value,value:row.children[2].value}));
  const related={}; if(!$('relatedPanel').hidden)relatedKeys.forEach(k=>{if($('related_'+k).value)related[k]=$('related_'+k).value;});
  const selected=Array.from($('columns').querySelectorAll('input:checked')).map(i=>i.value);
  if(!selected.length)throw Error('Select at least one column.');
  return {table:$('table').value,search:$('search').value,filters,related,columns:selected,sort,descending,limit:100};
}
async function apply(reset=true) {
  const id=++serial;
  $('status').className='';$('status').textContent='Searching…';$('export').disabled=true;$('saveQuery').disabled=true;$('previous').disabled=true;$('next').disabled=true;
  try {
    const spec=reset?specFromControls():{...active}; if(reset)offset=0;
    const response=await fetch(queryURL('query',{...spec,columns:[],offset}));
    const data=await response.json(); if(!response.ok)throw Error(data.error||'Query failed'); if(id!==serial)return;
    active=spec; total=data.total;
    render(data.rows,spec.columns);
    $('status').textContent=total.toLocaleString()+' matching rows';
    $('page').textContent=total?`${offset+1}–${Math.min(offset+data.rows.length,total)} of ${total.toLocaleString()}`:'No matches';
    $('previous').disabled=offset===0;$('next').disabled=offset+100>=total;
    $('export').disabled=false;$('saveQuery').disabled=false;
    history.replaceState(null,'','#spec='+encodeURIComponent(JSON.stringify(spec)));
  } catch(error) { if(id!==serial)return;$('status').className='error';$('status').textContent=error.message;$('previous').disabled=true;$('next').disabled=true; }
}
function drill(table, fields) {
  const valid=metadata.views[table].columns.map(c=>c.name);
  const filters=Object.entries(fields).filter(([k,v])=>valid.includes(k)&&v!==null&&v!==undefined).map(([column,value])=>({column,op:'eq',value}));
  configure(table,filters);apply();
}
function actions(row) {
  const cell=element('td'); cell.className='actions';
  const scope={sample_id:row.sample_id};
  if(row.orf_id){
    cell.append(button('ORF annotations',()=>drill('annotation_explorer',{...scope,orf_id:row.orf_id})),button('Pathway links',()=>drill('pathway_gene_explorer',{...scope,orf_id:row.orf_id})),button('ORF abundance',()=>drill('abundance_explorer',{...scope,feature_type:'orf',feature_id:row.orf_id})));
  }
  if(row.pathway_id)cell.append(button('Pathway genes',()=>drill('pathway_gene_explorer',{...scope,entity_id:row.entity_id,pathway_id:row.pathway_id})));
  if(row.contig_id)cell.append(button('Contig ORFs',()=>drill('orf_explorer',{...scope,contig_id:row.contig_id})));
  if(row.entity_id)cell.append(button('All mapped MAG ORFs',()=>drill('mag_orf_explorer',{...scope,entity_id:row.entity_id})),button('Entity pathways',()=>drill('pathway_explorer',{...scope,entity_id:row.entity_id})),button('MAG input genes',()=>drill('mag_gene_explorer',{...scope,entity_id:row.entity_id})));
  if(row.representative_orf_id)cell.append(button('Group members',()=>drill('orf_groups',{...scope,representative_orf_id:row.representative_orf_id})));
  if(row.annotation_id!==undefined)cell.append(button('EC / reaction terms',()=>drill('annotation_terms',{annotation_id:row.annotation_id})));
  if(row.source_id!==undefined)cell.append(button('Source record',()=>drill('sources',{source_id:row.source_id})));
  if(row.sample_id&&!row.orf_id&&!row.contig_id&&!row.pathway_id&&!row.entity_id)cell.append(button('Sample ORFs',()=>drill('orf_explorer',scope)));
  if(row.path||row.source){
    const path=row.path||row.source;
    if(typeof path==='string'&&!path.startsWith('/')&&!path.split('/').includes('..')) {
      const a=element('a','Open source file');a.href='../'+path.split('/').map(encodeURIComponent).join('/');a.target='_blank';a.rel='noopener';cell.append(a);
    }
  }
  return cell;
}
function render(rows,selected) {
  const head=$('results').querySelector('thead'),body=$('results').querySelector('tbody'); head.replaceChildren();body.replaceChildren();
  const hr=element('tr');hr.append(element('th','Related results'));
  selected.forEach(c=>{const th=element('th');th.append(button(c+(sort===c?(descending?' ↓':' ↑'):''),()=>{descending=sort===c?!descending:false;sort=c;apply();}));hr.append(th);});head.append(hr);
  rows.forEach(row=>{const tr=element('tr');tr.append(actions(row));selected.forEach(c=>tr.append(element('td',row[c]===null?'':String(row[c]))));body.append(tr);});
}
function saveJSON() {
  const blob=new Blob([JSON.stringify({schema_version:metadata.schema_version,report_generated_utc:metadata.generated_utc,query:active},null,2)],{type:'application/json'});
  const link=element('a');link.href=URL.createObjectURL(blob);link.download='metapathways-query.json';link.click();setTimeout(()=>URL.revokeObjectURL(link.href),1000);
}
function fromHash() {
  const params=new URLSearchParams(location.hash.slice(1));
  try {
    const spec=params.has('spec')?JSON.parse(params.get('spec')):null;
    configure(spec?.table||params.get('table')||'samples',spec?.filters||[]);
    if(spec){$('search').value=spec.search||'';relatedKeys.forEach(k=>$('related_'+k).value=spec.related?.[k]||'');sort=spec.sort||'';descending=!!spec.descending; if(spec.columns)$('columns').querySelectorAll('input').forEach(i=>i.checked=spec.columns.includes(i.value));}
    apply();
  } catch(error){$('status').textContent='Cannot read this saved query: '+error.message;}
}
async function init() {
  if(location.protocol==='file:'){$('offline').hidden=false;$('status').textContent='Start the local portal to search these results.';document.querySelector('.controls').hidden=true;return;}
  try {
    const response=await fetch(api('meta'));if(!response.ok)throw Error('Open the URL printed by metapathways report --serve.');metadata=await response.json();
    Object.entries(metadata.views).forEach(([key,value])=>$('table').append(option(key,`${value.label} (${value.rows.toLocaleString()})`)));
    const run=metadata.run_details||{};
    const details=[run.mp_version?'MP '+run.mp_version:'Report built with MP '+metadata.report_mp_version,
      metadata.sample_paths.length+' samples',run.command,run.executor,run.status].filter(Boolean);
    $('runDetails').textContent=details.join(' · ');
    $('generated').textContent='Report updated '+new Date(metadata.generated_utc).toLocaleString();
    $('table').onchange=()=>{configure($('table').value);apply();};$('addFilter').onclick=()=>filterRow();$('apply').onclick=()=>apply();$('reset').onclick=()=>{configure($('table').value);apply();};
    $('search').onkeydown=e=>{if(e.key==='Enter')apply();};$('previous').onclick=()=>{offset=Math.max(0,offset-100);apply(false);};$('next').onclick=()=>{offset+=100;apply(false);};
    $('export').onclick=()=>{const a=element('a');a.href=queryURL('export',active);a.download='metapathways-subset.csv';a.click();};$('saveQuery').onclick=saveJSON;
    window.addEventListener('hashchange',fromHash);fromHash();
  }catch(error){$('status').className='error';$('status').textContent=error.message;}
}
init();
