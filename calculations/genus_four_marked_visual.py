"""Self-contained exact-value explorer; no server or external assets."""
import json
from pathlib import Path
OUT=Path('results/genus-four-marked/2026-11-06')
data=json.loads((OUT/'explorer-data.json').read_text())
compact=[{'name':d['geometry'],'id':d['slug'],'kind':d['kind'],
          'rows':[[r['beta'],r['H_degree'],r['alpha_0'],r['comparison_sign'],r['twisted_unmarked'],r['ordinary_unmarked']] for r in d['rows']]}
         for d in data['examples']]
page=r'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<style>
:root{--font-sans:ui-sans-serif,system-ui,sans-serif;--font-mono:ui-monospace,monospace;--foreground:#e7e9ee;--muted-foreground:#a7afbd;--card:#20252e;--border:#46505d;--primary:#9eb8ff;--chart-1:#9eb8ff;--chart-2:#46505d;--radius:10px}
body{margin:0;color:var(--foreground);font:14px/1.5 var(--font-sans)}p{margin:0 0 12px}label{display:block;font-size:12px;color:var(--muted-foreground)}select{box-sizing:border-box;width:100%;padding:7px;border:1px solid var(--border);border-radius:6px;background:var(--card);color:var(--foreground);font:inherit;margin:4px 0 12px}option{background:var(--card)}.grid{display:grid;grid-template-columns:1fr 1fr;gap:0 12px}.card{background:var(--card);padding:16px;border-radius:var(--radius);margin:4px 0 18px}.values{display:grid;grid-template-columns:1fr 1fr;gap:12px}.value{font:16px var(--font-mono);overflow-wrap:anywhere;margin:5px 0 12px}.muted,small{color:var(--muted-foreground)}svg{display:block;width:100%;height:168px}svg text{fill:var(--foreground);font:12px var(--font-sans)}.legend{display:flex;gap:18px;font-size:12px}.dot{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:5px}h3{font-size:15px;margin:0 0 8px}@media(max-width:480px){.values{grid-template-columns:1fr}.value{font-size:15px}}
</style></head><body>
<p><strong>Genus 4 · up to two labeled markings</strong><br><span class="muted">43 presentations · degree bound 2 · exact fractions</span></p>
<label>Ambient presentation<select id="target"></select></label>
<div class="grid"><label>Curve class<select id="beta"></select></label><label>Insertions<select id="profile"></select></label></div>
<div class="card"><div class="values"><div><small>Ambient twisted coefficient Z</small><div class="value" id="twisted"></div></div><div><small>Ordinary compact GW number N</small><div class="value" id="ordinary"></div></div></div><div class="muted" id="formula"></div><div class="muted" id="aggregate"></div></div>
<h3>Independent relations found within the Pixton span</h3>
<label>Chow codimension<select id="degree"><option value="2">2</option><option value="3" selected>3</option></select></label>
<div class="legend"><span><i class="dot" style="background:var(--chart-1)"></i>Found in this search</span><span><i class="dot" style="background:var(--chart-2)"></i>Full generated Pixton span</span></div>
<svg id="chart" viewBox="0 0 650 168" role="img" aria-label="Relation ranks by number of markings"></svg>
<p class="muted" style="font-size:12px">Class search: codimensions 0–3 only. Top numerical coefficients are calculated separately in dimensions 9–11. An unknown compact value is displayed as unknown; numerical equality is not used as a Pixton certificate.</p>
<script>
const data=__DATA__,ranks=__RANKS__;
const el=id=>document.getElementById(id);const profiles=[[],['1'],['psi'],['H'],...['1','psi','H'].flatMap(a=>['1','psi','H'].map(b=>[a,b]))];
data.forEach((d,i)=>el('target').add(new Option(`${d.id} · ${d.name}`,i)));profiles.forEach((p,i)=>el('profile').add(new Option(p.length?p.map(x=>x==='psi'?'ψ':x).join(', '):'unmarked',i)));
function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){[a,b]=[b,a%b]}return a}
function frac(s){let p=String(s).split('/');return[BigInt(p[0]),BigInt(p[1]||1)]}
function scaled(s,m){if(!m)return'0';if(s===null)return'unknown';let[a,b]=frac(s);a*=BigInt(m);let c=gcd(a,b);return b/c===1n?String(a/c):`${a/c}/${b/c}`}
function difference(a,b){let[x,y]=frac(a),[u,v]=frac(b);let num=x*v-u*y,den=y*v,c=gcd(num,den);return den/c===1n?String(num/c):`${num/c}/${den/c}`}
function update(){let d=data[+el('target').value],row=d.rows[+el('beta').value],p=profiles[+el('profile').value];if(!row)return;let k=p.filter(x=>x==='psi').length,h=p.filter(x=>x==='H').length,m=p.includes('1')?0:row[1]**h;for(let i=0;i<k;i++)m*=6+h+i;let z=scaled(row[4],m),ordinary=scaled(row[5],m),power=row[2]+k+h-9-p.length;el('twisted').textContent=z;el('ordinary').textContent=ordinary;el('formula').textContent=`H·β = ${row[1]} · marking multiplier ${m} · top t power ${power}`;let b=power<0?scaled(z,-1):power===0&&ordinary!=='unknown'?difference(scaled(ordinary,+row[3]),z):'undetermined';el('aggregate').textContent=`CJR sign ${row[3]} · numerical infinity aggregate: ${b}`}
function changeTarget(){el('beta').replaceChildren();let d=data[+el('target').value];d.rows.forEach((r,i)=>el('beta').add(new Option('('+r[0].join(', ')+')',i)));update()}
function chart(){let r=+el('degree').value,rows=ranks.filter(x=>x.codimension===r),max=Math.max(...rows.map(x=>x.pixton_rank)),left=82,space=430;el('chart').innerHTML=rows.map((d,i)=>{let y=12+i*51;return`<text x="0" y="${y+19}">${d.markings} marks</text><rect x="${left}" y="${y}" width="${space*d.pixton_rank/max}" height="28" rx="3" fill="var(--chart-2)"/><rect x="${left}" y="${y+5}" width="${space*d.combined_rank/max}" height="18" rx="2" fill="var(--chart-1)"/><text x="${left+space*d.pixton_rank/max+10}" y="${y+19}">${d.combined_rank} / ${d.pixton_rank}</text>`}).join('')}
el('target').onchange=changeTarget;el('beta').onchange=update;el('profile').onchange=update;el('degree').onchange=chart;
el('target').value=data.findIndex(d=>d.id==='G15');changeTarget();el('beta').value=2;el('profile').value=11;update();chart();
</script></body></html>'''
page=page.replace('__DATA__',json.dumps(compact,separators=(',',':'))).replace('__RANKS__',json.dumps(data['relations']['degrees'],separators=(',',':')))
(OUT/'explorer.html').write_text(page)
print(OUT/'explorer.html')
