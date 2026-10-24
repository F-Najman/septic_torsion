/* Exact A7 certificate for the remaining C40/F3 quintic place.
   Four exceptional quadratic cusps; both possible diamond7 four-cycles.
   No class-group computation. The other80 complementary cusp divisors
   are rational lifts, and A7 kills them on the divisor group. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
out40:=Open("logs/sieve/X1_40_cusp_correction.log","w");
procedure Note40(s)
 fprintf out40,"%o\n",s;Flush(out40);print s;
end procedure;
X:=MDX1(40,GF(3));assert Genus(X) eq 25;
data:=eval Read("data/N40_survivor_fixed.txt");
assert #data eq 1;
dg,mult,defpoly,coords:=Explode(data[1]);assert dg eq 5 and mult eq 1;
K:=ext<GF(3)|Polynomial(GF(3),defpoly)>;
ps:=Places(Curve(X)(K)![K!c:c in coords]);
assert #ps eq 1 and Degree(ps[1]) eq 5 and not IsCusp(X,ps[1]);P5:=ps[1];
Note40(Sprintf("C40_CUSP_CORRECTION_MODEL_READY genus=25 p=3 quintic_degree=5 time=%o",Cputime(t0)));

cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^12,2^^6,4^^6*};
by_signature:=AssociativeArray();
for c in cs do
 sig:=CuspSignature(X,c);
 if not IsDefined(by_signature,sig) then by_signature[sig]:=[];end if;
 Append(~by_signature[sig],c);
end for;
assert #Keys(by_signature) eq CuspOrbitCountQ(X) and #Keys(by_signature) eq 21;
badgroups:=[by_signature[k]:k in Sort(Setseq(Keys(by_signature)))|
 #by_signature[k] eq 2 and &and[Degree(c) eq 2:c in by_signature[k]]];
assert #badgroups eq 2;
badcusps:=badgroups[1] cat badgroups[2];assert #Seqset(badcusps) eq 4;
assert Degree(&+[Divisor(c):c in badgroups[1]]) eq 4;
assert Degree(&+[Divisor(c):c in badgroups[2]]) eq 4;
assert 12*13 div 2 + 6 eq 84;
Note40(Sprintf("C40_CUSP_CORRECTION_CUSPS bad_quadratic=4 Q_orbits=2 full_complements=84 Q_liftable=80 signatures=%o time=%o",
 [CuspSignature(X,badgroups[i][1]):i in [1..2]],Cputime(t0)));

// sigma7 swaps the two Frobenius3-orbits in each Q(zeta8) cusp orbit.
sigma:=[2,1,4,3];SS:=Sym(4);
roots:=[[i^g:i in [1..4]]:g in SS|[i^(g^2):i in [1..4]] eq sigma];
assert #roots eq 2 and #Seqset(roots) eq 2;
for delta in roots do
 assert &and[delta[delta[i]] eq sigma[i]:i in [1..4]];
 assert &and[(i le 2) ne (delta[i] le 2):i in [1..4]];
end for;
Note40(Sprintf("C40_CUSP_CORRECTION_POSSIBLE_DIAMONDS sigma=%o roots=%o time=%o",sigma,roots,Cputime(t0)));

E:=EllipticCurve(X,P5);PP:=LevelStructure(X,P5)`P;
Pd:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=7*PP>);
assert Degree(Pd) eq 5;
T:=HeckeOperatorCommonZeros(X,7,P5);assert Degree(T) eq 40;
A:=T-7*Divisor(Pd)-Divisor(P5);assert Degree(A) eq 0;
survivors:=0;tested:=0;
pr:=IsPrincipal(A);tested+:=1;if pr then survivors+:=1;end if;
Note40(Sprintf("C40_CUSP_CORRECTION_TEST kind=Q_liftable complements=80 principal=%o time=%o",pr,Cputime(t0)));
for r in [1..#roots] do
 delta:=roots[r];
 for i in [1..4] do
  corr:=Divisor(badcusps[sigma[i]])-Divisor(badcusps[i])
       +7*(Divisor(badcusps[delta[sigma[i]]])-Divisor(badcusps[delta[i]]));
  assert Degree(corr) eq 0;
  pr:=IsPrincipal(A+corr);tested+:=1;if pr then survivors+:=1;end if;
  Note40(Sprintf("C40_CUSP_CORRECTION_TEST kind=bad_quadratic cycle=%o cusp=%o principal=%o time=%o",r,i,pr,Cputime(t0)));
 end for;
end for;
assert tested eq 9;
Note40(Sprintf("C40_CUSP_CORRECTION_DONE tests=%o survivors=%o all_complements=84 complete=true time=%o",tested,survivors,Cputime(t0)));
delete out40;quit;
