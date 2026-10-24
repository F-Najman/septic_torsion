/* Necessary rational class condition for all cusp complements, before
   constructing any characteristic-zero cusp or Riemann--Roch space.
   Requires rank zero and rational cuspidal generation separately proved. */
SetColumns(0);SetMemoryLimit(3*10^9);t0:=Cputime();
if Type(N) eq MonStgElt then N:=StringToInteger(N);end if;
if Type(p) eq MonStgElt then p:=StringToInteger(p);end if;
if Type(D) eq MonStgElt then D:=StringToInteger(D);end if;
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX1(N,GF(p));C:=Curve(X);cusps:=Cusps(X);
orbits:=[Type(O) eq PlcCrvElt select Divisor(O) else O:O in CuspOrbitsQ(X)];
c0:=Support([O:O in orbits|Degree(O) eq 1][1])[1];assert IsCusp(X,c0);
Cl,m1,m2:=ClassGroup(C);Free:=FreeAbelianGroup(#orbits);
gens:=[m2(O-Degree(O)*Divisor(c0)):O in orbits];
h:=hom<Free->Cl|gens>;G:=Image(h);
print "RATIONAL_CUSPIDAL_IMAGE",N,p,Invariants(G),Cputime(t0);
survivors:=[];
for l in Split(Read(SURV),"\n") do
 if #l eq 0 then continue;end if;
 data:=eval l;dv:=DivisorGroup(C)!0;
 for item in data do
  dg,mult,defpoly,coords:=Explode(item);
  K:=ext<GF(p)|Polynomial(GF(p),defpoly)>;
  pt:=C(K)![K!c:c in coords];pls:=Places(pt);assert #pls eq 1;
  assert Degree(pls[1]) eq dg;dv+:=mult*Divisor(pls[1]);
 end for;
 Append(~survivors,dv);
end for;
complements:=AssociativeArray();
for k in {D-Degree(s):s in survivors} do
 assert k ge 0;
 if k eq 0 then complements[k]:=[DivisorGroup(C)!0];continue;end if;
 result:=[];
 for lambda in Partitions(k) do
  partial:=[DivisorGroup(C)!0];
  for e in Sort(Setseq(Seqset(lambda))) do
   ce:=[Divisor(c):c in cusps|Degree(c) eq e];
   if #ce eq 0 then partial:=[];break;end if;
   mm:=#[t:t in lambda|t eq e];
   choices:=[&+[z:z in S]:S in Multisets(Seqset(ce),mm)];
   partial:=[a+b:a in partial,b in choices];
  end for;
  result cat:=partial;
 end for;
 complements[k]:=result;
 print "CUSP_COMPLEMENTS",k,#result;
end for;
total:=0;passing:=0;
for i in [1..#survivors] do
 dv:=survivors[i];cs:=complements[D-Degree(dv)];hits:=0;
 for cc in cs do
  total+:=1;
  if m2(dv+cc-D*Divisor(c0)) in G then
   hits+:=1;print "CLASS_FILTER_HIT",i,cc;
  end if;
 end for;
 passing+:=hits;
 print "CLASS_FILTER_SURVIVOR",i,"degree",Degree(dv),"complements",#cs,"hits",hits,Cputime(t0);
end for;
printf "DEGREE7_CLASS_FILTER_DONE N=%o p=%o D=%o inputs=%o tested=%o passing=%o complete=true\n",N,p,D,#survivors,total,passing;
quit;
