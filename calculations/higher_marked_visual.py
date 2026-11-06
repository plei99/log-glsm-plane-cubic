"""A self-contained chart and exact three-mark value explorer."""
import json
from pathlib import Path
OUT=Path('results/higher-marked/2026-11-06')
d=json.loads((OUT/'explorer-data.json').read_text())
page=r'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><style>
:root{--foreground:#e7e9ee;--muted-foreground:#a7afbd;--card:#20252e;--border:#46505d;--chart-1:#a6baff;--chart-2:#4c5665;--radius:10px;--font-sans:system-ui,sans-serif;--font-mono:ui-monospace,monospace}
body{margin:0;color:var(--foreground);font:14px/1.5 var(--font-sans)}p{margin:0 0 12px}h3{font-size:15px;margin:8px 0}label,small,.muted{color:var(--muted-foreground)}label{font-size:12px;display:block}select{box-sizing:border-box;width:100%;background:var(--card);color:var(--foreground);border:1px solid var(--border);border-radius:6px;padding:7px;margin:4px 0 12px;font:inherit}.grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}.card{background:var(--card);border-radius:var(--radius);padding:16px;margin:4px 0 16px}.value{font:16px var(--font-mono);overflow-wrap:anywhere;margin:5px 0 10px}svg{width:100%;height:205px;display:block}svg text{font:12px var(--font-sans);fill:var(--foreground)}.legend{font-size:12px}.dot{display:inline-block;width:9px;height:9px;margin-right:5px;border-radius:2px}@media(max-width:450px){.values{grid-template-columns:1fr}}
</style></head><body><p><strong>Genera 2 and 3 · up to three marked points</strong><br><span class="muted">43 presentations · degree ≤ 2 · exact fractions</span></p>
<div class="grid"><label>Genus<select id="genus"><option value="2">2</option><option value="3">3</option></select></label><label>Curve class<select id="beta"></select></label></div>
<label>Ambient presentation<select id="target"></select></label><label>Insertion profile<select id="profile"></select></label>
<div class="card"><div class="grid values"><div><small>Ambient twisted coefficient</small><div class="value" id="z"></div></div><div><small>Ordinary compact GW number</small><div class="value" id="N"></div></div></div><div class="muted" id="factor"></div></div>
<h3>Codimension-two relations as markings increase</h3><div class="legend"><span class="dot" style="background:var(--chart-1)"></span>Found here &nbsp; <span class="dot" style="background:var(--chart-2)"></span>Generated Pixton span</div><svg id="chart" viewBox="0 0 650 205" role="img" aria-label="Relation ranks at zero through three marked points"></svg>
<p class="muted" id="beyond"></p><p class="muted" style="font-size:12px">Chow search bounded to codimensions 0–2. Top numerical coefficients are computed separately. Unknown compact values remain unknown; all found relation spans have exact Pixton certificates.</p>
<script>
const data=__DATA__,summary=__SUMMARY__;const el=id=>document.getElementById(id);
const profiles=[[]];for(let n=1;n<=3;n++){let a=[[]];for(let i=0;i<n;i++)a=a.flatMap(x=>['1','psi','H'].map(y=>[...x,y]));profiles.push(...a)}
profiles.forEach((p,i)=>el('profile').add(new Option(p.length?p.map(x=>x==='psi'?'ψ':x).join(', '):'unmarked',i)));
let active=[];function gcd(a,b){a=a<0n?-a:a;while(b){[a,b]=[b,a%b]}return a}
function multiply(s,m){if(!m)return'0';if(s===null)return'unknown';let x=s.split('/'),a=BigInt(x[0])*BigInt(m),b=BigInt(x[1]||1),c=gcd(a,b);return b/c===1n?String(a/c):`${a/c}/${b/c}`}
function values(){let x=active[+el('target').value],r=x.rows[+el('beta').value],p=profiles[+el('profile').value],g=+el('genus').value,k=p.filter(x=>x==='psi').length,h=p.filter(x=>x==='H').length,m=p.includes('1')?0:r[1]**h;for(let j=0;j<k;j++)m*=2*g-2+h+j;el('z').textContent=multiply(r[4],m);el('N').textContent=multiply(r[5],m);el('factor').textContent=`H·β = ${r[1]} · multiplier ${m} · top t power ${r[2]+k+h-(3*g-3+p.length)}`}
function beta(){el('beta').replaceChildren();active[+el('target').value].rows.forEach((r,i)=>el('beta').add(new Option('('+r[0].join(', ')+')',i)));values()}
function genus(){let old=active[+el('target').value]?.id;active=data.filter(x=>x.genus===+el('genus').value);el('target').replaceChildren();active.forEach((x,i)=>el('target').add(new Option(`${x.id} · ${x.name}`,i)));let chosen=active.findIndex(x=>x.id===(old||'G15'));el('target').value=chosen<0?0:chosen;beta();chart()}
function chart(){let g=+el('genus').value,rows=summary.degrees.filter(x=>x.genus===g&&x.codimension===2),max=Math.max(...rows.map(x=>x.pixton_rank));el('chart').innerHTML=rows.map((x,i)=>{let y=12+48*i,l=75,w=430;return`<text x="0" y="${y+19}">${x.markings} marks</text><rect x="${l}" y="${y}" width="${w*x.pixton_rank/max}" height="28" rx="3" fill="var(--chart-2)"/><rect x="${l}" y="${y+5}" width="${w*x.found_rank/max}" height="18" rx="2" fill="var(--chart-1)"/><text x="${l+w*x.pixton_rank/max+10}" y="${y+19}">${x.found_rank} / ${x.pixton_rank}</text>`}).join('');let b=summary.unit_pullback_comparison.find(x=>x.genus===g&&x.markings===3&&x.codimension===2);el('beyond').textContent=`At three marks: ${b.found_directions_beyond_unit_pullbacks} found directions beyond unit pullbacks of the two-mark relations.`}
el('genus').onchange=genus;el('target').onchange=beta;el('beta').onchange=values;el('profile').onchange=values;genus();el('beta').value=2;el('profile').value=profiles.findIndex(p=>p.join(',')==='psi,psi,psi');values();
</script></body></html>'''
page=page.replace('__DATA__',json.dumps(d['examples'],separators=(',',':'))).replace('__SUMMARY__',json.dumps(d['relations'],separators=(',',':')))
(OUT/'explorer.html').write_text(page);print(OUT/'explorer.html')
