/* New degree7 branch. The rational T2 annihilator and the lower-degree
   single-place tests are inherited degree-independent inputs. */
SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("hecke_two.m");
X:=MDX1(49,GF(3));assert Genus(X) eq 69;
places:=NoncuspidalPlaces(X,7);
reps:=PlacesUpToDiamondFast(X,places);
print "N49_DEGREE7_SETUP",#places,#reps,Cputime(t0);
remaining:=0;
for i in [1..#reps] do
    print "N49_DEGREE7_CASE_START",i,#reps,Cputime(t0);
    x:=reps[i];T:=HeckeTwoCommonZeros(X,x);
    A:=T-2*DiamondOperator(X,2,Divisor(x))-Divisor(x);
    assert Degree(A) eq 0;
    principal:=IsPrincipal(A);
    print "N49_DEGREE7_CASE",i,"principal",principal,Cputime(t0);
    if principal then remaining+:=1;end if;
end for;
printf "N49_DEGREE7_DONE total=%o remaining=%o complete=true time=%o\n",#reps,remaining,Cputime(t0);
quit;
