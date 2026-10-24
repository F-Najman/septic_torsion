/* Degree55 pure good7 branch. Allbad2 is excluded separately by rankzero and gonality>14. Ordinary Hecke3 annihilates rational Jacobian torsion. No cusp enumeration is required here. */
SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("hecke_two.m");
X:=MDX1(55,GF(2));assert Genus(X) eq 81;
places:=NoncuspidalPlaces(X,7);
reps:=PlacesUpToDiamondFast(X,places);
print "N55_DEGREE7_SETUP",#places,#reps,Cputime(t0);
remaining:=0;
for i in [1..#reps] do
    print "N55_DEGREE7_CASE_START",i,#reps,Cputime(t0);
    x:=reps[i];T:=HeckeOperatorCommonZeros(X,3,Divisor(x));
    A:=T-3*DiamondOperator(X,3,Divisor(x))-Divisor(x);
    assert Degree(A) eq 0;
    principal:=IsPrincipal(A);
    print "N55_DEGREE7_CASE",i,"principal",principal,Cputime(t0);
    if principal then remaining+:=1;end if;
end for;
printf "N55_DEGREE7_DONE total=%o remaining=%o complete=true time=%o\n",#reps,remaining,Cputime(t0);
quit;
