"""Split an expensive target over canonical profiles, then merge exact spans."""
from log_glsm.core.hodge import _load_vendored_admcycles
import argparse,json,gzip,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from sage.all import *
from calculations import higher_marked_classes as common
def merge(index,g):
    t=common.genus_four_targets()[index];folder=common.OUT/'parts'/(t.slug+'-g'+str(g));records=[]
    profiles=common.higher_marked_profiles(canonical=True)
    for i,profile in enumerate(profiles):
        with gzip.open(folder/('part%s.json.gz'%i),'rt') as f:data=json.load(f)
        assert data['profile_subset']==[list(profile)] and data['genus']==g
        records.append(data)
    states={(n,r):{'echelon':{},'rows':[],'origins':[],'candidates':0} for n in range(4) for r in range(3)}
    adm=_load_vendored_admcycles()
    for part,data in enumerate(records):
        for item in data['relations']:
            n,r=item['markings'],item['codimension'];s=states[n,r];s['candidates']+=item['candidates']
            columns=len(adm.TautologicalRing(g,n).generators(r))
            for row,origin in zip(item['rows'],item['origins']):
                v=vector(QQ,columns,{int(j):QQ(c) for j,c in row.items()});rem=copy(v)
                for pivot,b in sorted(s['echelon'].items()):
                    if rem[pivot]:rem-=rem[pivot]*b
                if not rem:continue
                pivot=next(j for j,x in enumerate(rem) if x);s['echelon'][pivot]=rem/rem[pivot]
                s['rows'].append(row);s['origins'].append(dict(origin,split_part=part))
    result=dict(records[0],profile_subset=None,split_part=None,split_over_profiles=len(profiles),
        finite_classes=[row for data in records for row in data['finite_classes']],
        seconds=sum(data['seconds'] for data in records),seconds_measure='sum of per-profile worker wall times',
        relations=[{'markings':n,'codimension':r,'candidates':s['candidates'],'rows':s['rows'],'origins':s['origins']} for (n,r),s in states.items()])
    path=common.OUT/(t.slug+'-g'+str(g)+'-classes.json.gz');temp=path.with_suffix('.tmp')
    with gzip.open(temp,'wt') as f:json.dump(result,f)
    temp.replace(path);print('MERGED',path,[(n,r,len(s['rows'])) for (n,r),s in states.items() if s['rows']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,default=42);p.add_argument('--genus',type=int,default=2)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--parts',type=int,nargs='*');p.add_argument('--merge-only',action='store_true');args=p.parse_args()
    if args.merge_only:merge(args.index,args.genus)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(common.scan,args.index,args.genus,2,(profile,),i)
                     for i,profile in enumerate(common.higher_marked_profiles(canonical=True))
                     if args.parts is None or i in args.parts]
            for f in as_completed(futures):print('PART DONE',f.result(),flush=True)
        print('All profile parts complete; merge with --merge-only after other writers have stopped.',flush=True)
