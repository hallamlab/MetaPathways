'use strict';
const $=id=>document.getElementById(id);
function el(tag,text){const node=document.createElement(tag);if(text!==undefined)node.textContent=text;return node;}
function url(kind,spec){const u=new URL('../api/'+kind,location.href);if(spec)u.searchParams.set('spec',JSON.stringify(spec));return u.href;}
let metadata, reviewed=null;
function spec(){return {mode:document.querySelector('input[name="mode"]:checked').value,selections:[...document.querySelectorAll('#selections input:checked')].filter(c=>!c.disabled).map(c=>({sample_id:c.dataset.sample,...(c.dataset.mag?{entity_id:c.dataset.mag}:{})}))};}
function changed(){reviewed=null;$('download').disabled=true;$('review').textContent='';const n=spec().selections.length;$('status').textContent=n+' selections';$('preview').disabled=!n;}
function checkbox(sample,mag,label){const line=el('label'),box=el('input');box.type='checkbox';box.dataset.sample=sample;if(mag)box.dataset.mag=mag;line.append(box,document.createTextNode(' '+label));line.style.display='block';return [line,box];}
function fileSize(bytes){
  if(bytes<1024)return bytes+' B';
  const units=['KiB','MiB','GiB','TiB'];let value=bytes/1024,index=0;
  while(value>=1024&&index<units.length-1){value/=1024;index++;}
  return value.toFixed(1)+' '+units[index];
}
function reviewText(data){
  const lines=[(data.mode==='tables'?'EDA tables only':'Complete dataset')+' — '+data.total_files+' files in the ZIP',''];
  data.items.forEach(item=>{
    const c=item.file_counts;
    lines.push(item.scope.sample_id+(item.scope.entity_id?' / '+item.scope.entity_id:' / whole sample')+' — '+item.total_files+' files');
    lines.push('  '+c.table_csv+' CSV tables (including provenance)');
    if(c.sequence_fasta)lines.push('  '+c.sequence_fasta+' generated sequence FASTA files');
    if(c.annotation_subset)lines.push('  '+c.annotation_subset+' subset GFF/GenBank files');
    if(data.mode==='complete')lines.push('  '+c.retained_original+' copied source files ('+fileSize(item.retained_bytes)+' uncompressed; copied files only)');
    lines.push('');
  });
  lines.push(data.metadata_files+' shared metadata files: manifest.json and README.txt');
  return lines.join('\n');
}
async function init(){
  if(location.protocol==='file:'){$('status').textContent='Start metapathways report -o OUTPUT --serve --no-rebuild to use downloads.';return;}
  try{
    const response=await fetch(url('meta'));if(!response.ok)throw Error('Open the local URL printed by metapathways report --serve.');metadata=await response.json();
    for(const sample of metadata.scope_samples){
      const group=el('details');group.open=true;group.append(el('summary',sample));
      const [wholeLabel,whole]=checkbox(sample,null,'Whole sample (includes every MAG)');group.append(wholeLabel);
      const mags=el('div');mags.style.marginLeft='1.5rem';const boxes=[];
      const all=el('button','Select all MAGs in this sample');all.type='button';group.append(all);
      for(const mag of metadata.scope_mags.filter(m=>m.sample_id===sample)){
        const [label,box]=checkbox(sample,mag.entity_id,mag.entity_id+(mag.has_contigs?'':' — contig membership unavailable'));
        box.dataset.available=mag.has_contigs?'yes':'no';box.disabled=!mag.has_contigs;box.onchange=changed;boxes.push(box);mags.append(label);
      }
      whole.onchange=()=>{boxes.forEach(b=>{b.disabled=whole.checked||b.dataset.available==='no';if(whole.checked)b.checked=false;});changed();};
      all.onclick=()=>{whole.checked=false;whole.onchange();boxes.forEach(b=>{if(!b.disabled)b.checked=true;});changed();};
      group.append(mags);$('selections').append(group);
    }
    $('clear').onclick=()=>{document.querySelectorAll('#selections input').forEach(b=>{b.checked=false;b.disabled=b.dataset.available==='no';});changed();};
    $('allSamples').onclick=()=>{$('clear').click();document.querySelectorAll('#selections input:not([data-mag])').forEach(b=>{b.checked=true;b.onchange();});changed();};
    $('allMags').onclick=()=>{$('clear').click();document.querySelectorAll('#selections input[data-mag]').forEach(b=>{if(!b.disabled)b.checked=true;});changed();};
    document.querySelectorAll('input[name="mode"]').forEach(b=>b.onchange=changed);
    $('preview').onclick=async()=>{const current=spec(),key=JSON.stringify(current);$('status').textContent='Reviewing selection…';try{const r=await fetch(url('download-preview',current));const data=await r.json();if(!r.ok)throw Error(data.error);if(key!==JSON.stringify(spec()))return;reviewed=current;$('review').textContent=reviewText(data);$('status').textContent=data.note;$('download').disabled=false;}catch(e){$('status').textContent=e.message;}};
    $('download').onclick=()=>{if(!reviewed)return;const a=el('a');a.href=url('download',reviewed);a.download='metapathways-selected-results.zip';a.target='_blank';a.rel='noopener';a.click();$('status').textContent='Archive requested. Your browser will begin downloading when preparation finishes.';};changed();
  }catch(e){$('status').textContent=e.message;}
}
init();
