#!/usr/bin/env python3
"""Audit complete exact winding logs and degree-seven local inequalities."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,re,ast
H=Path(__file__).resolve().parent
TARGETS=[81,99,115,161,253,299,391,437,529]
def ds(n):return [d for d in range(1,n+1) if n%d==0]
def phi(n):return sum(math.gcd(a,n)==1 for a in range(1,n+1))
def order(n,a):
 k=1;t=a%n
 while t not in [1,n-1]:t=t*a%n;k+=1
 return k
rows=[];pending=[]
ff=(H/'logs/finite_fields/torsion_F2_degree7.log').read_text()+(H/'logs/finite_fields/torsion_F2_99.log').read_text()+(H/'logs/finite_fields/torsion_F2_81.log').read_text()
for n in TARGETS:
 logs=[H/f'logs/winding/X1_{n}.log']
 logs=[p for p in logs if p.exists() and 'WINDING_CHARACTER_PROBE_DONE' in p.read_text()]
 if not logs:pending.append(n);continue
 p=logs[0];s=p.read_text()
 assert not any(x in s for x in ['Runtime error','System error','Assertion failed','User error'])
 lines=re.findall(r'WINDING N=(\d+) order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+) time=([0-9.]+)',s)
 assert lines
 if n==529:
  replay=json.loads((H/'data/winding529_krylov_replay.json').read_text())
  assert replay['verified'] and replay['Krylov_rank']==replay['dimension']==45 and replay['splitting_prime']==2531
  assert 'WINDING_GOOD_PRIME p=2531 all_denominators_invertible=true' in s
  assert 'WINDING_KRYLOV N=529 order=253 span_rank=45 dimension=45 verified=true' in s
 chars=[dict(zip(['N','order','conductor','dimension','winding_dimension','coefficient_degree'],map(int,t[:6])),CPU_seconds=float(t[6])) for t in lines]
 assert all(t['N']==n and (t['order']==1 or t['winding_dimension']==t['dimension']) for t in chars)
 fin=re.findall(r'WINDING_FINAL N=(\d+) genus=(\d+) total_dimension=(\d+) bad_orders=(\[.*?\]) admissible=(\[.*?\])',s)
 assert len(fin)==1
 nn,g,dd=map(int,fin[0][:3]);bad=ast.literal_eval(fin[0][3]);admiss=ast.literal_eval(fin[0][4]);assert nn==n and g==dd and bad in [[],[1]] and 2 in admiss
 assert sum(t['dimension']*phi(t['order']) for t in chars)==g
 I=Fraction(n*n,2)
 for q in ds(n):
  if q>1 and len(ds(q))==2:I*=Fraction(q*q-1,q*q)
 assert I.denominator==1
 gg=1+I/12-Fraction(sum(phi(d)*phi(n//d) for d in ds(n)),4)
 assert gg==g
 G=Fraction(325*I,32768);o=order(n,2)
 assert G>28 and 7%o
 assert re.search(rf'TORSION n={n} p=2 : k with Z/{n} over F_p\^k = \[\s*\]',ff)
 rows.append({'N':n,'a':2,'diamond_order':o,'index':int(I),'genus':g,'gonality_lower_bound':str(G),'finite_field_good_degrees_through7':[], 'only_deficient_character_orders':bad,'all_nontrivial_characters_full_winding':True,'character_orbits':chars,'log':p.name,'log_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
out={'degree':7,'excluded_cyclic_orders':[r['N'] for r in rows],'pending_exact_winding_logs':pending,'rows':rows,'method':'dyadic finite image diamond; CM septic curves excluded by Bourdon-Pollack Theorem1.2','status':'all listed numerical, winding and finite-field inputs verified; proof is in new_large_exclusions.tex'}
(H/'data/new_large_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('NEW_LARGE_CERTIFICATE_DONE',out['excluded_cyclic_orders'],'pending',pending)
for r in rows:print(r['N'],r['diamond_order'],r['index'],r['genus'],r['gonality_lower_bound'])
