#!/usr/bin/env python3
"""Exact degree7 numerical/rank-provenance checks; no BSD assumptions.

Newform records are existing proved LMFDB inputs from the quintic paper.
This recomputes the complete dimension sum at every divisor level and
records every nontrivial character requiring a separate exact kernel test.
"""
from fractions import Fraction
import hashlib,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
INPUTS=HERE/'data'
cache=json.loads((INPUTS/'lmfdb_newforms_cache.json').read_text())
rows=[(119,3,7,96),(121,2,8,110),(125,3,7,50),(133,64,8,36),
      (143,67,8,120),(169,3,7,52),(187,122,8,160),(209,34,8,180),
      (221,3,7,96),(247,144,8,72),(289,3,7,136),(323,101,8,288),
      (361,3,7,114)]
local=[(77,3,60),(85,2,32),(91,3,24),(95,11,72)]
def divs(n):return [d for d in range(1,n+1) if n%d==0]
def phi(n):return sum(math.gcd(a,n)==1 for a in range(1,n+1))
def primes(n):return [p for p in divs(n) if p>1 and len(divs(p))==2]
def index(n):
 x=Fraction(n*n,2)
 for p in primes(n):x*=Fraction(p*p-1,p*p)
 assert x.denominator==1
 return int(x)
def order(n,a):
 v=a%n;t=1
 while v not in (1,n-1):v=v*a%n;t+=1
 return t
allrows=[];chars=[]
for item in rows+[(n,a,None,cm) for n,a,cm in local]:
 n,a,k,cm=item;idx=index(n)
 genus=1+idx//12-sum(phi(d)*phi(n//d) for d in divs(n))//4
 forms=[];dim=0
 for m in divs(n):
  assert str(m) in cache
  mult=len(divs(n//m))
  for f in cache[str(m)]:
   assert f['analytic_rank_proved'] is True and f['analytic_rank'] is not None
   dim+=f['dim']*mult
   if f['analytic_rank']>0:forms.append({'level':m,'multiplicity':mult,**f})
 assert dim==genus
 nontriv=sorted({(f['level'],f['char_orbit_label']) for f in forms if f['char_order']>1})
 assert all(m==n for m,lab in nontriv)
 lower=Fraction(325*idx,32768)
 assert 7%order(n,a)!=0 and cm>7
 if k is None:assert lower>28
 else:assert n%3 and lower>7*k
 allrows.append({'N':n,'a':a,'method':'local2' if k is None else 'global3','k3':k,
                 'index':idx,'gonality_lower_bound':str(lower),'diamond_order':order(n,a),
                 'CM_minimum':cm,'genus':genus,'newform_dimension_sum':dim,
                 'positive_rank_forms':forms,'nontrivial_positive_characters':nontriv})
 chars.append((n,a,[lab for m,lab in nontriv]))
text=['SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);',
      'AttachSpec("mdmagma/Magma/magma.spec");']
for n,a,labs in chars:
 text.append(f'N:={n};a:={a};')
 for lab in labs:
  text.append('idx:=[j:j in [1..N-1]|GCD(j,N) eq 1 and ConreyCharacterOrbitLabel(N,j) eq Sprintf("%o.%o",N,"'+lab+'")];')
  text.append('assert #idx gt 0 and &and[ConreyCharacterValue(N,j,a) eq 1:j in idx];')
  text.append('printf "SEPTIC_CHARACTER_KERNEL N=%o a=%o orbit=%o indices=%o verified=true\\n",N,a,"'+lab+'",idx;')
 text.append('printf "SEPTIC_FINITE_IMAGE_INPUT N=%o a=%o nontrivial_positive_characters=%o verified=true\\n",N,a,'+str(len(labs))+';')
text.extend(['print "SEPTIC_INHERITED_CHARACTER_CERTIFICATE_DONE";','quit;'])
(HERE/'check_inherited_characters.m').write_text('\n'.join(text)+'\n')
result={'degree':7,'rank_input':'data/lmfdb_newforms_cache.json',
        'rank_input_sha256':hashlib.sha256((INPUTS/'lmfdb_newforms_cache.json').read_bytes()).hexdigest(),
        'excluded_cyclic_orders':sorted(r['N'] for r in allrows),'rows':allrows,
        'status':'numerical_and_rank_inputs_verified; separate_character_and_finite_field_logs_required'}
(HERE/'data/inherited_large_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print('SEPTIC_LARGE_NUMERICAL_DONE',len(allrows),result['excluded_cyclic_orders'])
