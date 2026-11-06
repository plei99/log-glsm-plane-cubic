"""Inline explorer for the bounded toric calculation and exact span comparison."""
import json
from pathlib import Path
out=Path('results/toric-complete-intersections/2026-11-06')
data={'catalog':json.loads((out/'catalog.json').read_text()),
      'union':json.loads((out/'union-summary.json').read_text()),
      'ordinary':json.loads((out/'ordinary-cy3-gw.json').read_text())}
html='''<!doctype html><html><head><meta charset="utf-8"><style>
*{box-sizing:border-box}body{margin:0;color:var(--foreground);font:14px var(--font-sans)}p{line-height:1.5;margin:8px 0 14px}.muted{color:var(--muted-foreground)}
.controls{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}label{display:flex;flex-direction:column;gap:6px;font-size:12px}
select{font:inherit;color:var(--foreground);background:var(--card);border:1px solid var(--border);border-radius:var(--radius);padding:7px;max-width:100%}
.formula{font:18px var(--font-mono);margin:14px 0}.values{display:flex;flex-wrap:wrap;gap:12px;margin:14px 0}.metric{background:var(--card);padding:16px;border-radius:var(--radius);flex:1;min-width:145px}.metric .number{font:19px var(--font-mono);margin-top:9px;overflow-wrap:anywhere}
.chart{height:205px;display:flex;flex-direction:column;gap:17px;margin:18px 0}.row{display:grid;grid-template-columns:40px 1fr 100px;gap:12px;align-items:center}.track{height:16px;background:var(--secondary);border-radius:4px;overflow:hidden}.fill{height:16px;background:var(--chart-1)}.value{text-align:right;font-variant-numeric:tabular-nums}.chips{display:flex;gap:6px;flex-wrap:wrap}.chip{padding:8px 11px;border-radius:var(--radius);border:1px solid var(--border);font-size:12px}.negative{color:var(--chart-1);border-color:var(--chart-1)}.constant{color:var(--chart-2);border-color:var(--chart-2)}
</style></head><body>
<p><strong>Toric complete-intersection calculations</strong><br><span class="muted">Twelve targets, genera two and three, and Mori tuples with coefficient sum ≤ 2. S denotes the negative section of F₁; F is a fiber.</span></p>
<div class="controls"><label>Target<select id="target"></select></label><label>Genus<select id="genus"><option>2</option><option selected>3</option></select></label><label>Mori tuple β<select id="beta"></select></label></div>
<p id="description"></p><div class="formula" id="grading"></div><div class="chips" id="chips"></div>
<div class="values"><div class="metric"><span class="muted">Ambient twisted ∫Z</span><div class="number" id="twisted"></div></div><div class="metric"><span class="muted">Ordinary CY3 N</span><div class="number" id="ordinary"></div></div><div class="metric"><span class="muted">Infinity aggregate ∫B</span><div class="number" id="infinity"></div></div></div>
<p class="muted" id="explanation"></p>
<p><strong>Genus-three relation span</strong><br><span class="muted">Blue: toric rank. Full bar: earlier projective-family rank in the same codimension. Exact certificates place the entire toric span inside the earlier span.</span></p>
<div class="chart" id="chart"></div><p class="muted">The toric union has 19 genus-two and 187 genus-three relations. Adding these to the earlier search leaves its combined rank at 463: zero additional independent directions. Every relation is certified inside Pixton’s system.</p>
<script>const DATA=__DATA__;
const names={T1:'(1,1) in P¹ × P²',T2:'(2,3) in P¹ × P²',T3:'(1,2) in P¹ × P³',T4:'(2,4) in P¹ × P³',T5:'(1,2) in P² × P²',T6:'(3,3) in P² × P²',T7:'(1,1),(1,1) in P² × P²',T8:'(1,1),(2,2) in P² × P²',T9:'S+F+H in F₁ × P¹',T10:'2S+3F+2H in F₁ × P¹',T11:'S+F+H in F₁ × P²',T12:'2S+3F+3H in F₁ × P²'};
const select=document.getElementById('target');DATA.catalog.forEach((s,i)=>{const o=document.createElement('option');o.value=i;o.textContent=s.target.id+' · '+names[s.target.id];select.appendChild(o);});select.value=11;
function neg(s){return s==='0'?'0':s.startsWith('-')?s.slice(1):'-'+s;}
function choose(){const s=DATA.catalog[+select.value];const b=document.getElementById('beta');b.replaceChildren();Object.keys(s.genera[0].top_twisted_integrals).forEach(k=>{const o=document.createElement('option');o.value=k;o.textContent='('+k+')';b.appendChild(o);});if(s.target.id==='T12')b.value='2,0,0';update();}
function update(){const s=DATA.catalog[+select.value],t=s.target,g=+document.getElementById('genus').value,key=document.getElementById('beta').value,beta=key.split(',').map(Number),gs=s.genera.find(x=>x.genus===g);let degree=0,k=0;t.ambient.factors.forEach(f=>{if(f==='F1'){degree+=(t.c1[k+1]-t.c1[k])*beta[k]+t.c1[k]*beta[k+1];k+=2;}else{degree+=t.c1[k]*beta[k];k++;}});const alpha=t.dimension*(g-1)-degree,dim=3*g-3,power=alpha-dim,Z=gs.top_twisted_integrals[key];
document.getElementById('description').textContent=t.kind+' '+(t.dimension===2?'surface':'threefold')+' · Euler characteristic '+t.chern_numbers.euler+' · '+(t.ambient.factors[0]==='F1'?'β = (dS, dF, dH)':'β = (d₁, d₂)');document.getElementById('grading').textContent='α = '+t.dimension+' × '+(g-1)+' − '+degree+' = '+alpha;
const chips=document.getElementById('chips');chips.replaceChildren();for(let r=0;r<=dim;r++){const p=alpha-r,d=document.createElement('div');d.className='chip '+(p<0?'negative':p===0?'constant':'muted');d.innerHTML='r='+r+' · t<sup>'+p+'</sup>';chips.appendChild(d);}
const known=DATA.ordinary.find(x=>x.target===t.id&&x.genus===g&&x.beta.join(',')===key);document.getElementById('twisted').textContent=Z;document.getElementById('ordinary').textContent=t.kind==='CY3'?(known?known.ordinary_gw:'unknown'):'see compact class';
let B='unspecified';if(power<0||power===0&&t.kind==='K3')B=neg(Z);else if(power===0&&known)B=known.infinity_aggregate_integral;document.getElementById('infinity').textContent=B;
document.getElementById('explanation').textContent='The values are top integrals on M̄'+g+'. Blue powers are negative-t CJR identities; the cyan t⁰ term is the signed compact comparison. '+(t.kind==='CY3'?'CY3s have no negative-t range here. Their ordinary GW values are supplied only where separately established.':t.kind==='K3'?'Ordinary K3 theory vanishes, so the t⁰ aggregate also equals −Z.':'Ordinary Fano classes at t⁰ remain explicit in the report.');}
select.addEventListener('input',choose);document.getElementById('genus').addEventListener('input',update);document.getElementById('beta').addEventListener('input',update);choose();
const chart=document.getElementById('chart');DATA.union.degrees.filter(x=>x.genus===3&&x.previous_rank).forEach(x=>{const d=document.createElement('div');d.className='row';d.innerHTML='<span>r='+x.codimension+'</span><div class="track"><div class="fill" style="width:'+(100*x.rank/x.previous_rank)+'%"></div></div><span class="value">'+x.rank+' / '+x.previous_rank+'</span>';chart.appendChild(d);});
</script></body></html>'''.replace('__DATA__',json.dumps(data,separators=(',',':')))
(out/'toric-explorer.html').write_text(html)
