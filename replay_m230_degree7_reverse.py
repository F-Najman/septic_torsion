#!/usr/bin/env python3
"""Independent integer-lattice replay of the septimic A3-first sieve."""
import json,math,re
from collections import Counter
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
HERE=Path(__file__).resolve().parent

def read_expressions(path):
    s=path.read_text().replace('[*','[').replace('*]',']')
    dec=json.JSONDecoder();out=[];pos=0
    while pos<len(s):
        while pos<len(s) and s[pos].isspace():pos+=1
        if pos==len(s):break
        val,end=dec.raw_decode(s,pos);out.append(val);pos=end
    return out

ar=read_expressions(HERE/'data/M230_D7_reverse_arithmetic.txt');assert len(ar)==1
inv,relations,hgens,degrees,cusps,records,local,survivors=ar[0]
assert inv==[4,4,4,136,16320,0] and len(records)==368 and len(local)==64 and survivors==[]
allvectors=relations+hgens+cusps+[v for r in records for v in r[2]]+[v for r in local for v in r[4]]
assert all(len(v)==6 and v[-1]==0 for v in allvectors)
B=hermite_normal_form(Matrix([r[:-1] for r in relations+hgens]).T)
assert B.shape==(5,5) and abs(int(B.det()))==4352
assert all(B[i,j]==0 for i in range(5) for j in range(i))
b=[[int(B[i,j]) for j in range(5)] for i in range(5)]

def key(v):
    v=list(v[:5]);assert len(v)==5
    for i in range(4,-1,-1):
        q=v[i]//b[i][i]
        for j in range(i+1):v[j]-=q*b[j][i]
    assert all(0<=v[i]<b[i][i] for i in range(5))
    return tuple(v)
zero=(0,)*5
assert all(key(v)==zero for v in relations+hgens)
cc=[key(v) for v in cusps];C=[{zero}]
for k in range(1,6):
    C.append({key([x+y for x,y in zip(v,c)])
              for d,c in zip(degrees,cc) if d<=k for v in C[k-d]})
assert list(map(len,C))==[1,9,41,67,68,68]
allowed=[{key([-a for a in v]) for v in s} for s in C]
checks=0
for idx,e,sums,passed,parts in records:
    assert 1<=idx<=1404 and 2<=e<=7 and passed is False and sums
    assert sum(d*m for d,m,n in parts)==e
    assert len(sums)<=math.prod(n for d,m,n in parts)
    for v in sums:
        checks+=1;assert key(v) not in allowed[7-e]
for idx,d,split,image_degrees,classes in local:
    assert 1<=idx<=64 and 2<=d<=7
    if split:
        assert len(classes)==len(image_degrees)==3
        assert all(d%f==0 for f in image_degrees)
    else:
        assert len(classes)==len(image_degrees)==1 and 3*d<=7
        assert (3*d)%image_degrees[0]==0
rrs=read_expressions(HERE/'data/A3_D7_RR_records_v2.txt')
# Probe and final sieve each wrote one complete expression to this append-only diagnostic.
assert len(rrs) in (1,2) and all(r==rrs[0] for r in rrs)  # the committed file holds two identical copies
rr=rrs[-1];assert [r[0] for r in rr]==list(range(1,257))
az=read_expressions(HERE/'data/A3_D7_group_arithmetic.txt');assert len(az)==1
zinv,zrel,zhs,zdegrees=az[0]
assert zinv==[2,2,2,2,2,2,102,10608,0]
assert sorted(zdegrees)==[1,1,1,1,2,2,4,4]
assert all(len(v)==9 and v[-1]==0 for v in zrel+zhs+[r[1] for r in rr])
Z=hermite_normal_form(Matrix([v[:-1] for v in zrel+zhs]).T)
ZR=hermite_normal_form(Matrix([v[:-1] for v in zrel]).T)
assert Z.shape==ZR.shape==(8,8)
assert abs(int(ZR.det()))==69249024 and abs(int(Z.det()))==270504
assert abs(int(ZR.det()))//abs(int(Z.det()))==256

def zkey(v,basis):
    v=list(v[:8])
    for i in range(7,-1,-1):
        q=v[i]//int(basis[i,i])
        for j in range(i+1):v[j]-=q*int(basis[j,i])
    return tuple(v)
assert all(zkey(r[1],Z)==(0,)*8 for r in rr)
assert len({zkey(r[1],ZR) for r in rr})==256

hist=Counter();total=nc=0
for idx,h,dim,patterns in rr:
    hist[dim]+=1;assert len(patterns)==(0 if dim==0 else (7**dim-1)//6)
    for ds,ms,noncusp in patterns:
        assert len(ds)==len(ms) and all(d>0 for d in ds) and all(m>0 for m in ms)
        assert sum(d*m for d,m in zip(ds,ms))==7
        total+=1;nc+=noncusp
assert hist=={0:84,1:108,2:48,3:16} and total==1404 and nc==816
log=(HERE/'logs/m2/X1_2_30_filter.log').read_text()
assert not re.search('error|interrupt|segmentation|failed',log,re.I)
assert 'M230_D7_REVERSE_DONE target_classes=256 target_divisors=1404 noncuspidal=816 inert_rejected=448 Yclass_rejected=368 remaining=0 complete=true' in log
result={'status':'passed','classes':256,'target_divisors':1404,'noncuspidal_targets':816,
        'inert_rejected':448,'class_rejected':368,'explicit_class_sum_checks':checks,
        'cyclic_quotient_index':4352,'A3_image_order':256,'A3_quotient_index':270504,'cusp_complement_cardinalities':list(map(len,C)),
        'remaining':0}
(HERE/'data/M230_D7_reverse_arithmetic.json').write_text(json.dumps(ar[0],indent=2)+'\n')
(HERE/'data/A3_D7_RR_records_v2.json').write_text(json.dumps(rr,indent=2)+'\n')
(HERE/'data/M230_D7_reverse_replay.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
