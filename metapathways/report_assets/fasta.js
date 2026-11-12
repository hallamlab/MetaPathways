'use strict';
const $=id=>document.getElementById(id);
function url(kind,spec){const u=new URL('../api/'+kind,location.href);u.searchParams.set('spec',JSON.stringify(spec));return u.href;}
async function request(kind,spec){const r=await fetch(url(kind,spec));const data=await r.json();if(!r.ok)throw Error(data.error||'Sequence request failed');return data;}
async function init(){
  if(location.protocol==='file:'){$('status').textContent='Open this viewer through metapathways report --serve.';return;}
  try{
    const lookup=JSON.parse(new URLSearchParams(location.hash.slice(1)).get('lookup'));
    $('target').textContent=lookup.sample_id+' / '+(lookup.feature_id||lookup.contig_id);
    const data=await request('sequence-options',lookup);let serial=0;
    for(const record of data.records){const option=document.createElement('option');option.value=record.sequence_record_id;option.textContent=record.sequence_type+' · '+record.sequence_id+' · '+record.length.toLocaleString()+' residues · '+record.source;$('record').append(option);}
    $('record').disabled=!data.records.length;
    $('status').textContent=data.records.length?'Choose one available sequence to inspect.'+(data.more_records?' More than 500 source records exist; use Downloads for the complete set.':''):'No retained sequence could be linked to this feature or contig. See Import notes and Downloads.';
    $('record').onchange=async()=>{const id=++serial;$('fasta').textContent='';$('download').hidden=true;$('source').textContent='';if(!$('record').value)return;try{
      const record=data.records.find(r=>String(r.sequence_record_id)===$('record').value);
      const spec={sample_id:lookup.sample_id,sequence_record_id:record.sequence_record_id};$('status').textContent='Loading sequence…';
      const result=await request('sequence-preview',spec);if(id!==serial)return;
      $('source').textContent='Source: '+record.source;$('fasta').textContent=result.fasta;
      $('status').textContent=result.truncated?'Preview: first 10,000 of '+result.length.toLocaleString()+' residues. Download the full FASTA to inspect the remainder.':'Complete sequence · '+result.length.toLocaleString()+' residues.';
      $('download').href=url('fasta',spec);$('download').download='metapathways-sequence.fasta';$('download').hidden=false;
    }catch(e){if(id===serial)$('status').textContent=e.message;}};
  }catch(e){$('status').textContent=e.message;}
}
init();
