/* New degree7 branch. The rational T2 annihilator and the lower-degree
   single-place tests are inherited degree-independent inputs. */
SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("hecke_two.m");
X:=MDX1(11,GF(3));
ps:=NoncuspidalPlaces(X,3);assert #ps ge 2;
for i in [1..2] do
 a:=HeckeTwoCommonZeros(X,ps[i]);b:=HeckeOperator(X,2,ps[i]);
 assert a eq b;print "HECKE2_INDEPENDENT_CHECK",i,"exact_divisors_equal=true";
end for;
print "HECKE2_INDEPENDENT_CHECK_DONE";quit;
