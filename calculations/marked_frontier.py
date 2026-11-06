"""Incremental higher-marking search, with exact orbit closure and RSS records."""
import argparse,fcntl,gc,gzip,json,os,resource,threading,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from itertools import combinations_with_replacement
from sage.all import *
from calculations import higher_marked_classes as common
from functools import lru_cache
from log_glsm.targets.convex_gkm import (
    convex_gkm_targets,
)

OUT=Path(os.environ.get('MARKED_FRONTIER_OUT','results/marked-frontier/2026-11-06'))


@lru_cache(maxsize=1)
def targets():return common.genus_four_targets()+convex_gkm_targets()


def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    with (gzip.open(temp,'wt') if path.suffix=='.gz' else temp.open('w')) as f:
        json.dump(data,f,default=lambda x:int(x) if x in ZZ else str(x))
    temp.replace(path)


def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)


class RelationSpan:
    """Maintain an invariant span using adjacent transpositions, not n! actions."""
    def __init__(self,g,n,r):
        self.g,self.n,self.r=g,n,r
        self.echelon={};self.rows=[];self.origins=[];self.candidates=0
        self.generators=[]
        for i in range(n-1):
            a=list(range(1,n+1));a[i],a[i+1]=a[i+1],a[i];self.generators.append(tuple(a))

    def add(self,v,origin):
        self.candidates+=1;pending=[(v,origin)]
        while pending:
            w,why=pending.pop();rem=copy(w)
            for pivot,row in sorted(self.echelon.items()):
                if rem[pivot]:rem-=rem[pivot]*row
            if not rem:continue
            pivot=next(i for i,c in enumerate(rem) if c)
            self.echelon[pivot]=rem/rem[pivot];self.rows.append(common.sparse(w));self.origins.append(why)
            for mapping in self.generators:
                pending.append((common.permute_vector(w,self.g,self.n,self.r,mapping),
                                dict(why,adjacent_transposition=mapping)))

    def record(self):
        return {'markings':self.n,'codimension':self.r,'rows':self.rows,
                'origins':self.origins,'candidate_seeds':self.candidates}


def profiles(n):return tuple(combinations_with_replacement(('1','psi','H'),n))


def direct_profiles(n,include_units=False):return tuple(p for p in profiles(n) if include_units or n<4 or '1' not in p)


def result_path(index,g,n,r,part=None):
    name=targets()[index].slug+'-g%s'%g
    return OUT/('n%s-r%s'%(n,r))/(name+('.json.gz' if part is None else '-part%s.json.gz'%part))


def scan(index,g,n,max_r=2,part=None,parts=1,include_units=False):
    start=time.time();target=targets()[index];path=result_path(index,g,n,max_r,part)
    # Do not overlap the largest graph catalogues at adjacent marking levels.
    # A lower-mark checkpoint also supplies the subsequent unit reconstruction.
    if index==42 and g==2 and n>=5:
        while not result_path(index,g,n-1,max_r).exists():time.sleep(5)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.with_suffix('.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        selected=direct_profiles(n,include_units) if part is None else direct_profiles(n,include_units)[part::parts]
        if path.exists():
            old=read(path)
            recorded=old.get('direct_profiles',old['profiles'])
            assert (recorded==[list(p) for p in selected] or (part is None and recorded==[list(p) for p in profiles(n)])) and old['max_codimension']==max_r
            return str(path),'reused'
        states={r:RelationSpan(g,n,r) for r in range(max_r+1)};reference={};finite=[]
        for base in (4,5):
            b=common.higher_marked_backend(target,g,base)
            for beta in common.genus_four_degrees(target):
                for profile in selected:
                    co=b.marked_coefficients(beta,profile,max_r)
                    for r,powers in co.items():
                        state=states[r];origin={'beta':list(map(int,beta)),'profile':profile,'base':base}
                        for power,terms in powers.items():
                            if power<0:state.add(common.vector_of(terms,g,n,r),dict(origin,kind='pole',epsilon=power))
                        v=common.vector_of(powers.get(0,{}),g,n,r)
                        if base==4:
                            reference[beta,profile,r]=v
                            if v:finite.append({'beta':list(map(int,beta)),'profile':profile,'codimension':r,'vector':common.sparse(v)})
                            # Equal entries in a canonical profile occur in blocks.
                            # Adjacent swaps within each block generate its stabilizer.
                            for mapping in state.generators:
                                if all(profile[i]==profile[j-1] for i,j in enumerate(mapping)):
                                    state.add(common.permute_vector(v,g,n,r,mapping)-v,
                                              dict(origin,kind='profile_symmetry',permutation=mapping))
                        else:state.add(v-reference[beta,profile,r],dict(origin,kind='ray_difference'))
                print(target.slug,g,n,'part',part,'base',base,'beta',beta,'ranks',
                      [len(s.rows) for s in states.values()],flush=True)
        data={'target':target.record(),'genus':g,'markings':n,'max_codimension':max_r,'curve_bound':2,
              'rays':[4,5],'profiles':selected,'part':part,'parts':parts,'finite_canonical_classes':finite,
              'relations':[s.record() for s in states.values()],
              'seconds':time.time()-start,'peak_worker_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'ordered_recovery':'Apply exact label permutations to canonical profiles; no factorial division.'}
        dump(path,data)
    return str(path),data['seconds'],data['peak_worker_rss_kib']


def merge(index,g,n,r,parts,include_units=False):
    records=[read(result_path(index,g,n,r,i)) for i in range(parts)]
    assert sorted(tuple(p) for d in records for p in d['profiles'])==sorted(direct_profiles(n,include_units))
    states={k:RelationSpan(g,n,k) for k in range(r+1)}
    for data in records:
        for item in data['relations']:
            k=item['codimension'];size=len(common.converter(g).adm.TautologicalRing(g,n).generators(k))
            for row,origin in zip(item['rows'],item['origins']):
                states[k].add(vector(QQ,size,{int(i):QQ(c) for i,c in row.items()}),origin)
    result=dict(records[0],profiles=direct_profiles(n,include_units),part=None,parts=parts,
                relations=[s.record() for s in states.values()],
                finite_canonical_classes=[x for d in records for x in d['finite_canonical_classes']],
                seconds=sum(d['seconds'] for d in records),seconds_measure='sum of worker wall times',
                peak_worker_rss_kib=max(d['peak_worker_rss_kib'] for d in records))
    dump(result_path(index,g,n,r),result)


def monitor(stop,path,root=None):
    """RSS sum includes shared pages; MemAvailable and PSI measure system impact."""
    start=time.time();root=os.getpid() if root is None else root;path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w') as f:
        while True:
            table={}
            for p in Path('/proc').glob('[0-9]*'):
                try:
                    status={line.split(':')[0]:line.split(':')[1].strip() for line in (p/'status').read_text().splitlines() if ':' in line}
                    table[int(p.name)]=(int(status['PPid']),int(status.get('VmRSS','0 kB').split()[0]))
                except (OSError,KeyError,ValueError):pass
            descendants={root};changed=True
            while changed:
                extra={pid for pid,(parent,_) in table.items() if parent in descendants}-descendants
                changed=bool(extra);descendants.update(extra)
            mem={k:int(v.split()[0]) for k,v in (line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())}
            vm=dict(line.split() for line in Path('/proc/vmstat').read_text().splitlines())
            row={'seconds':time.time()-start,'task_rss_sum_kib':sum(table.get(p,(0,0))[1] for p in descendants),
                 'processes':len(descendants),'available_kib':mem['MemAvailable'],'total_kib':mem['MemTotal'],
                 'swap_used_kib':mem['SwapTotal']-mem['SwapFree'],
                 'swapin_pages':int(vm['pswpin']),'swapout_pages':int(vm['pswpout']),
                 'oom_kill':int(vm.get('oom_kill',0)),
                 'memory_pressure':Path('/proc/pressure/memory').read_text() if Path('/proc/pressure/memory').exists() else None}
            f.write(json.dumps(row)+'\n');f.flush()
            if root not in table or stop.wait(5):break


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--markings',type=int,required=True)
    p.add_argument('--codimension',type=int,default=2);p.add_argument('--workers',type=int,default=6)
    p.add_argument('--indices',type=int,nargs='*',default=list(range(54)))
    p.add_argument('--genera',type=int,nargs='*',default=[2,3]);p.add_argument('--direct-only',action='store_true')
    p.add_argument('--include-unit-profiles',action='store_true',help='Localize units directly as well, for full-profile replay.')
    args=p.parse_args()
    assert 0<=args.markings<=8 and 0<=args.codimension<=3
    jobs=[];splits=[]
    for i in args.indices:
        for g in args.genera:
            path=result_path(i,g,args.markings,args.codimension)
            if path.exists():
                data=read(path)
                if args.include_unit_profiles and data.get('direct_profiles',data['profiles'])!=[list(p) for p in profiles(args.markings)]:
                    raise ValueError('Full-profile replay requires a fresh MARKED_FRONTIER_OUT; this checkpoint used unit reconstruction.')
                continue
            count=5 if i==42 and g==2 else 1
            if count>1:splits.append((i,g,args.markings,args.codimension,count,args.include_unit_profiles))
            jobs += [(i,g,args.markings,args.codimension,k if count>1 else None,count,args.include_unit_profiles) for k in range(count)]
    stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,OUT/('memory-n%s-%s.jsonl'%(args.markings,int(time.time())))))
    thread.start()
    try:
        with ProcessPoolExecutor(max_workers=args.workers,max_tasks_per_child=1) as pool:
            for f in as_completed([pool.submit(scan,*j) for j in jobs]):print('DONE',f.result(),flush=True)
        for spec in splits:merge(*spec)
        if args.markings>=4 and not args.direct_only:
            from calculations.frontier_forgetting import complete
            for i in args.indices:
                for g in args.genera:complete(i,g,args.markings,args.codimension)
    finally:stop.set();thread.join()
