"""Render the exact scan's grading explorer and relation-rank comparison."""
import json
from pathlib import Path
out=Path('results/complete-intersections/2026-11-06')
catalog=json.loads((out/'catalog.json').read_text())
union=json.loads((out/'union-summary.json').read_text())
data={'targets':catalog,'union':union}
html='''<!doctype html><html><head><meta charset="utf-8"><style>
*{box-sizing:border-box}body{margin:0;color:var(--foreground);font:14px var(--font-sans)}
p{line-height:1.5;margin:8px 0 14px}.muted{color:var(--muted-foreground)}
.controls{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}label{display:flex;flex-direction:column;gap:6px;font-size:12px}
select,input{font:inherit;color:var(--foreground);background:var(--card);border:1px solid var(--border);border-radius:var(--radius);padding:7px}
select{max-width:100%}.equation{font-size:20px;font-family:var(--font-mono);margin:16px 0 12px}
.chips{display:flex;flex-wrap:wrap;gap:7px;margin:12px 0}.chip{padding:9px 12px;border-radius:var(--radius);border:1px solid var(--border);min-width:80px}
.chip .small{font-size:11px;margin-top:5px}.negative{border-color:var(--chart-1);color:var(--chart-1)}
.constant{border-color:var(--chart-2);color:var(--chart-2)}.positive{color:var(--muted-foreground)}
.chart{height:250px;margin:18px 0 4px;display:flex;flex-direction:column;gap:15px}
.row{display:grid;grid-template-columns:40px 1fr 110px;gap:12px;align-items:center}.track{height:16px;background:var(--secondary);border-radius:4px;overflow:hidden}.fill{height:16px;background:var(--chart-1)}
.value{text-align:right;font-variant-numeric:tabular-nums}.legend{display:flex;gap:18px;flex-wrap:wrap;font-size:12px;margin:10px 0 18px}.dot{display:inline-block;width:8px;height:8px;border-radius:50%;margin-right:5px}
</style></head><body>
<p><strong>Which CJR coefficient gives a relation?</strong><br><span class="muted">Choose a target, genus and ambient curve degree. α = m(g − 1) − ID; the codimension-r term has power t<sup>α−r</sup>.</span></p>
<div class="controls"><label>Target<select id="target"></select></label><label>Genus<select id="genus"><option>2</option><option selected>3</option></select></label><label>Curve degree D<input id="degree" type="number" min="0" max="30" value="3" style="width:90px"></label></div>
<div class="equation" id="equation"></div><div class="chips" id="chips"></div>
<div class="legend"><span><i class="dot" style="background:var(--chart-1)"></i>Negative power: B = −Z</span><span><i class="dot" style="background:var(--chart-2)"></i>t⁰: compact GW comparison</span><span class="muted">Positive power: no automatic vanishing</span></div>
<p class="muted" id="targetnote"></p>
<p><strong>Genus-three relation space found across all 14 targets</strong><br><span class="muted">Independent strata expressions, with exact Pixton certificates. Each bar is scaled to the generated Pixton rank in that codimension.</span></p>
<div class="chart" id="chart"></div><p class="muted" id="totals"></p>
<script>const DATA=__DATA__;
const select=document.getElementById('target');
DATA.targets.forEach((s,i)=>{const o=document.createElement('option');o.value=i;o.textContent=s.target.name.replace('(5,)','(5)');select.appendChild(o);if(s.target.slug==='p5-5')select.value=i;});
function update(){const t=DATA.targets[+select.value].target;const g=+document.getElementById('genus').value;const D=Math.max(0,Math.min(30,Math.trunc(+document.getElementById('degree').value||0)));const alpha=t.dimension*(g-1)-t.index*D;const dim=3*g-3;
document.getElementById('equation').textContent='α = '+t.dimension+' × '+(g-1)+' − '+t.index+' × '+D+' = '+alpha;
const chips=document.getElementById('chips');chips.replaceChildren();for(let r=0;r<=dim;r++){const p=alpha-r;const div=document.createElement('div');div.className='chip '+(p<0?'negative':p===0?'constant':'positive');div.innerHTML='r = '+r+'<div class="small">t<sup>'+p+'</sup></div>';chips.appendChild(div);}
document.getElementById('targetnote').textContent=(t.index===0?'For this CY3, α = 3g − 3 at every curve degree: there are no negative-t slots. Auxiliary torus poles still supply tautological relations. ':'The infinity term B is the sum of all contributions containing infinity data. Its t⁰ comparison retains the ordinary compact GW class. ')+'The grading applies beyond the computed range; this scan evaluates D = 0, 1, 2, 3.';}
for(const id of ['target','genus','degree'])document.getElementById(id).addEventListener('input',update);update();
const chart=document.getElementById('chart');DATA.union.degrees.filter(x=>x.genus===3&&x.pixton_rank>0).forEach(x=>{const row=document.createElement('div');row.className='row';row.innerHTML='<span>r = '+x.codimension+'</span><div class="track"><div class="fill" style="width:'+(100*x.rank/x.pixton_rank)+'%"></div></div><span class="value">'+x.rank+' / '+x.pixton_rank+'</span>';chart.appendChild(row);});
const ranks=g=>DATA.union.degrees.filter(x=>x.genus===g).reduce((a,x)=>a+x.rank,0);document.getElementById('totals').textContent=ranks(2)+' independent genus-two relations and '+ranks(3)+' genus-three relations after removing repetitions across targets. Search: two auxiliary weight rays, degrees 0–3. All retained relations have exact rational Pixton certificates.';
</script></body></html>'''.replace('__DATA__',json.dumps(data,separators=(',',':')))
(out/'grading-and-relations.html').write_text(html)
