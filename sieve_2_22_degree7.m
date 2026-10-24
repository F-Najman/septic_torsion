/* New degree7 sieve; historical degree6 exclusions are not used.
   Rank0: DS Theorem4.1. gon_Q>=8: NV2602.03513v4 Table1.
   Odd-prime all-cusp criterion therefore applies withp3,d7.
   The only nonempty noncuspidal parts have single degree5,6,7 support.
   Exact cusp stabilizer data is degree-independent infrastructure. */
SetColumns(0);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("fast_diamond_hash_m2.m");
load "data/cusp_conditions_m2.m";
X:=MDX11(2,22,GF(3));C:=Curve(X);assert Genus(X) eq 16;
assert X`_E[1] eq 0 and X`_E[3] eq 0;
assert Sort([Degree(c):c in Cusps(X)]) eq TateCuspDegrees[<22,3>];
for k in [1,2] do
 assert 5 in TateCuspAllowed[<22,3>][k] and 7 in TateCuspAllowed[<22,3>][k];
end for;
places:=[NoncuspidalPlaces(X,d):d in [1..7]];
assert [#xs:xs in places] eq [0,0,0,0,60,100,240];
print "SEPTIC_M2_22_SETUP",Genus(X),[#xs:xs in places],Cputime(t0);
survivors:=[* *];tested:=0;
for d in [5..7] do
 rs:=PlacesUpToDiamondM2HashDirect(X,places[d]);assert #rs eq [12,20,48][d-4];
 print "SEPTIC_M2_22_REPS",d,#rs,Cputime(t0);passed:=0;
 for i in [1..#rs] do
  x:=rs[i];killed:=false;
  for q in [5,7] do
   T:=HeckeOperatorFast(X,q,x);delta:=DiamondOperator(X,q,x);
   A:=T-q*Divisor(delta)-Divisor(x);assert Degree(A) eq 0;
   principal:=IsPrincipal(A);
   print "SEPTIC_M2_22_CASE",d,i,q,principal,Cputime(t0);
   if not principal then killed:=true;break;end if;
  end for;
  tested+:=1;
  if not killed then
   passed+:=1;point:=RepresentativePoint(x);L:=GF(3^d);
   coords:=MDM2ToCommonField(Eltseq(point),L);
   def:=[Integers()!a:a in Coefficients(DefiningPolynomial(L,GF(3)))];
   values:=[[Integers()!a:a in Eltseq(z,GF(3))]:z in coords];
   Append(~survivors,[*d,i,def,values*]);
   print "SEPTIC_M2_22_SURVIVOR",d,i,def,values,Cputime(t0);
  end if;
 end for;
 print "SEPTIC_M2_22_TYPE",d,#rs,passed,Cputime(t0);
end for;
assert tested eq 80;
Write("data/sieve_2_22_degree7_survivors.m",
 "SepticM2N22Survivors := " cat Sprint(survivors) cat ";" : Overwrite:=true);
printf "SEPTIC_M2_22_DONE tested=%o remaining=%o complete=true time=%o\n",tested,#survivors,Cputime(t0);
quit;
