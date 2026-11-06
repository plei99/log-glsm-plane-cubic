"""Summarize host and guest samples without attributing all host activity to this task."""
import json,re,gzip
from pathlib import Path
OUT=Path('results/marked-frontier/2026-11-06')
baseline=json.loads((OUT/'host-baseline.json').read_text())
host=[json.loads(s) for s in (OUT/'host-memory.jsonl').read_text().splitlines()]
guest=[json.loads(s) for p in OUT.glob('memory-*.jsonl') for s in p.read_text().splitlines()]
workers=[]
for p in OUT.glob('n*-r*/*.json.gz'):
    with gzip.open(p,'rt') as f:d=json.load(f)
    workers.append({'file':str(p),'peak_rss_kib':d['peak_worker_rss_kib']})
peak=max(workers,key=lambda x:x['peak_rss_kib'])
result={'physical_gib':baseline['physical_bytes']/2**30,'guest_total_gib':max(x['total_kib'] for x in guest)/2**20,
        'minimum_guest_available_gib':min(x['available_kib'] for x in guest)/2**20,
        'largest_sampled_pool_rss_sum_gib':max(x['task_rss_sum_kib'] for x in guest)/2**20,
        'peak_recorded_worker':peak,'peak_worker_gib':peak['peak_rss_kib']/2**20,
        'guest_swapout_page_change':max(x['swapout_pages'] for x in guest)-min(x['swapout_pages'] for x in guest),
        'guest_swapin_page_change':max(x['swapin_pages'] for x in guest)-min(x['swapin_pages'] for x in guest),
        'guest_oom_kill_change':max(x['oom_kill'] for x in guest)-min(x['oom_kill'] for x in guest),
        'initial_compressor_gib':baseline['compressor_pages']*baseline['page_size']/2**30,
        'peak_host_compressor_gib':max(x['vm_stat']['Pages occupied by compressor']*x['page_size'] for x in host)/2**30,
        'host_swapout_gib_since_initial_snapshot':(host[-1]['vm_stat']['Swapouts']-baseline['swapout_pages'])*baseline['page_size']/2**30,
        'host_swapin_gib_since_initial_snapshot':(host[-1]['vm_stat']['Swapins']-baseline['swapin_pages'])*baseline['page_size']/2**30,
        'minimum_host_memory_free_percent':min(int(re.search(r'percentage: (\d+)',x['pressure']).group(1)) for x in host),
        'guest_psi_available':any(x['memory_pressure'] is not None for x in guest),
        'limitations':['RSS sums count shared pages more than once and are per pool, not a combined physical-memory total.',
                      'Host counters include other applications; they are not an exclusive measurement of this calculation.',
                      'The first guest monitor stopped because this kernel lacks PSI; later monitoring records that absence and continues.',
                      'Interrupted worker pools caused by a corrected Sage integer/decorator import error were not OOM events.']}
(OUT/'memory-summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
