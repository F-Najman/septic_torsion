/* Complete new degree7 rational-image sieve, independent of incomplete
   Hecke prefixes. N26/N28 have NC minimum degree5, so one support only.
   Every cusp complement has degree<=2; exact Tate Galois permutations
   show every such closed cusp divisor lifts overQ. */
SetColumns(0);SetMemoryLimit(8*10^9);t0:=Cputime();
if not assigned N then N:=26;elif Type(N) eq MonStgElt then N:=StringToInteger(N);end if;
assert N in [26,28];
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
Attach("fast_diamond_hash_m2.m");
load "data/cusp_conditions_m2_degree7.m";
X:=MDX11(2,N,GF(3));C:=Curve(X);assert Genus(X) eq 25;
cs:=Cusps(X);assert Sort([Degree(c):c in cs]) eq TateCuspDegrees[<N,3>];
for k in [1,2] do assert #TateCuspAllowed[<N,3>][k] eq EulerPhi(N);end for;
rat:= [c:c in cs|Degree(c) eq 1];c0:=rat[1];
places:=[NoncuspidalPlaces(X,d):d in [1..7]];
expected:=N eq 26 select [0,0,0,0,36,122,288] else [0,0,0,0,48,64,312];
assert [#ss:ss in places] eq expected;
reps:=[PlacesUpToDiamondM2HashDirect(X,places[d]):d in [5..7]];
er:=N eq 26 select [6,26,48] else [8,16,52];assert [#rs:rs in reps] eq er;
print "M2_D7_FULL_CLASS_SETUP",N,expected,er,Cputime(t0);
blocks:=[Divisor(c):c in cs|Degree(c) le 2];
if N eq 28 then
 FF:=FunctionField(C);b,c:=Explode([FF!f:f in X`_coordinates]);dd:=b*(1-b-c);
 EE:=EllipticCurve([0,c,0,dd,0]);QQ:=EE![b,b];fs:=[b,dd,c^2-4*dd];
 for k in [1..N div 2] do
  A:=k*QQ;assert A[3] ne 0 and A[1] ne 0;Append(~fs,A[1]);
 end for;
 groups:=AssociativeArray();
 for cusp in cs do
  sig:=[Valuation(f,cusp):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],cusp);else groups[sig]:=[cusp];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));
 unitsblocks:=[&+[Divisor(cusp):cusp in groups[sig]]:sig in sigs];
 blocks cat:=unitsblocks;
 print "M2_D7_FULL_CLASS_UNIT_BLOCKS",N,#sigs,[Degree(o):o in unitsblocks],Cputime(t0);
end if;
print "M2_D7_FULL_CLASS_GROUP_START",N,Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
J:=JH(2*N,[N+1]);assert Dimension(J) eq 25;
assert &*[a:a in Invariants(G)|a ne 0] eq Integers()!Evaluate(FrobeniusPolynomial(J,3),1);
hgens:=[m2(o-Degree(o)*Divisor(c0)):o in blocks];H:=sub<G|hgens>;
expectedH:=N eq 26 select [2,14,266,3990,11970] else [2,4,4,4,8,24,936,936];
assert Invariants(H) eq expectedH;
print "M2_D7_FULL_CLASS_GROUP",N,Invariants(G),Invariants(H),Cputime(t0);
records:=[* *];remaining:=[* *];tested:=0;
for pos in [1..3] do
 d:=pos+4;passed:=0;
 for i in [1..#reps[pos]] do
  x:=reps[pos][i];cl:=m2(Divisor(x)-d*Divisor(c0));inside:=cl in H;tested+:=1;
  Append(~records,[*d,i,Eltseq(cl),inside*]);
  if inside then passed+:=1;Append(~remaining,[*d,i,x*]);end if;
  print "M2_D7_FULL_CLASS_TEST",N,d,i,inside,Eltseq(cl),Cputime(t0);
 end for;
 print "M2_D7_FULL_CLASS_TYPE",N,d,#reps[pos],passed,Cputime(t0);
end for;
assert tested eq &+er;
export:=[* *];
for rec in remaining do
 d,i,x:=Explode(rec);pt:=RepresentativePoint(x);L:=GF(3^d);
 coords:=MDM2ToCommonField(Eltseq(pt),L);
 def:=[Integers()!a:a in Coefficients(DefiningPolynomial(L,GF(3)))];
 values:=[[Integers()!a:a in Eltseq(z,GF(3))]:z in coords];
 Append(~export,[*d,i,def,values*]);
end for;
Write(Sprintf("data/N%o_D7_full_class_arithmetic.txt",N),Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(a):a in hgens],records*]) : Overwrite:=true);
Write(Sprintf("data/N%o_D7_full_class_survivors.m",N),"M2D7FullClassSurvivors := " cat Sprint(export) cat ";" : Overwrite:=true);
printf "M2_D7_FULL_CLASS_DONE N=%o tested=%o remaining=%o complete=true time=%o\n",N,tested,#remaining,Cputime(t0);quit;
