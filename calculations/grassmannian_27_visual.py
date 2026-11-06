"""Display the geometrically justified line/conic contributions in degree two."""
import json
from fractions import Fraction
from pathlib import Path

out=Path('results/grassmannian-27/2026-11-06')
records=json.loads((out/'ordinary-gw.json').read_text())
data={}
for g,conic,line in ((0,Fraction(1225),Fraction(196,8)),
                     (2,Fraction(1225,240),Fraction(392,240)),
                     (3,Fraction(1225,6048),Fraction(1568,6048))):
    total=conic+line
    if g:assert str(total)==next(r['ordinary_gw'] for r in records if r['genus']==g and r['degree']==2)
    data[g]={'conics':str(conic),'lines':str(line),'total':str(total),'fraction':float(conic/total)}
html=r'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><style>
body{margin:0;font:14px/1.55 var(--font-sans,system-ui);color:var(--foreground)}*{box-sizing:border-box}.panel{background:var(--card);border-radius:var(--radius,12px);padding:16px}h3{margin:0 0 8px;font-size:16px}p{margin:10px 0}.muted{color:var(--muted-foreground)}select{font:inherit;color:var(--foreground);background:var(--secondary);border:1px solid var(--border);padding:6px 10px;border-radius:7px}.bar{display:flex;height:38px;border-radius:6px;overflow:hidden;margin:18px 0 10px;background:var(--secondary)}.conics{background:var(--chart-1)}.lines{background:var(--chart-2)}.columns{display:grid;grid-template-columns:1fr 1fr;gap:14px}.value{font:20px var(--font-mono,monospace);margin-top:4px}.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:6px}.total{margin-top:16px;padding-top:12px;border-top:1px solid var(--border)}.small{font-size:12px}
</style></head><body><div class="panel"><h3>Degree-two maps to Gr(2,7)[1⁷]</h3><p class="muted">A general target contains 196 isolated lines and 1225 isolated smooth conics. Their local contributions determine the full invariant in degree two.</p><label>Domain genus <select id="genus" onchange="draw()"><option value="3">3</option><option value="2">2</option><option value="0">0</option></select></label><div class="bar" role="img" id="bar"><div class="conics" id="conicbar"></div><div class="lines" id="linebar"></div></div><div class="columns"><div><span class="dot conics"></span>Degree-one maps onto conics<div class="value" id="conics"></div></div><div><span class="dot lines"></span>Degree-two maps onto lines<div class="value" id="lines"></div></div></div><div class="total">Ordinary connected invariant N<sub>g,2</sub> <span class="value" id="total"></span></div><p class="small muted">These stable-map contributions include contracted components and reducible domains. The two colors distinguish image curves, not domain topologies.</p></div><script>const data=__DATA__;function draw(){const g=document.getElementById('genus').value,d=data[g];for(const key of ['conics','lines','total'])document.getElementById(key).textContent=d[key];document.getElementById('conicbar').style.width=100*d.fraction+'%';document.getElementById('linebar').style.width=100*(1-d.fraction)+'%';document.getElementById('bar').setAttribute('aria-label',`Genus ${g}: conics contribute ${d.conics}; lines contribute ${d.lines}`)}draw();</script></body></html>'''
(out/'degree-two-contributions.html').write_text(html.replace('__DATA__',json.dumps(data)))
