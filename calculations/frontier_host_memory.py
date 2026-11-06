"""macOS host samples; run with system Python and stop after the task finishes."""
import json,re,subprocess,time
from pathlib import Path
out=Path(__file__).resolve().parents[1]/'results/marked-frontier/2026-11-06/host-memory.jsonl'
with out.open('a') as f:
    while True:
        vm=subprocess.check_output(['vm_stat'],text=True)
        page_size=int(re.search(r'page size of (\d+)',vm).group(1))
        stats={m.group(1):int(m.group(2)) for m in re.finditer(r'^([^:\n]+):\s+(\d+)\.',vm,re.M)}
        swap=subprocess.check_output(['sysctl','vm.swapusage'],text=True).strip()
        pressure=subprocess.check_output(['memory_pressure','-Q'],text=True).strip()
        f.write(json.dumps({'epoch':time.time(),'page_size':page_size,'vm_stat':stats,'swap':swap,'pressure':pressure})+'\n');f.flush()
        time.sleep(15)
