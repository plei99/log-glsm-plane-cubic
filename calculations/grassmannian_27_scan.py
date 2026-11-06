"""Run the Gr(2,7)[1^7] extension without changing the earlier sample outputs."""
from log_glsm.targets.grassmannian import grassmannian_seven_hyperplane_cy3
import shutil
from pathlib import Path
from types import SimpleNamespace
from sage.all import *
from calculations import grassmannian_ci_scan as scan
OUT=Path('results/grassmannian-27/2026-11-06')

if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    scan.OUT=OUT;scan.common.OUT=OUT;scan.common.ARGS=SimpleNamespace(max_degree=2)
    previous=Path('results/grassmannian-complete-intersections/2026-11-06')
    for name in ('stabilization-vectors.sobj','stabilization-cache.json'):
        if (previous/name).exists() and not (OUT/name).exists():shutil.copy2(previous/name,OUT/name)
    scan.common.prepare()
    target=grassmannian_seven_hyperplane_cy3()
    result=scan.scan_target(target.record())
    scan.dump(OUT/'catalog.json',[result])
    print('COMPLETE Gr(2,7)[1^7]',flush=True)
