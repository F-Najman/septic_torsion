SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_diamond_hash_m2.m");
Attach("fast_diamond_hash_cyclic_bc.m");
X:=MDX1(36,GF(5));xs:=NoncuspidalPlaces(X,3);assert #xs eq 48;
rs:=PlacesUpToDiamondCyclicHashBC(X,xs);legacy:=PlacesUpToDiamondFast(X,xs);
assert #rs eq #legacy;
K:=Parent(Evaluate(X`_coordinates[1],xs[1]));
allkeys:={MDM2FrobeniusKey([Evaluate(f,x):f in X`_coordinates],K,5,3):x in xs};
assert #allkeys eq #xs;seen:={[]};Exclude(~seen,[]);
for x in legacy do
 ks:=Seqset(MDCyclicBCHashKeys(X,x,K));
 actual:={MDM2FrobeniusKey([Evaluate(f,DiamondOperator(X,a,x)):f in X`_coordinates],K,5,3):a in [1..18]|GCD(a,36) eq 1};
 assert ks eq actual and IsEmpty(seen meet ks);seen join:=ks;
end for;
assert seen eq allkeys;
printf "CYCLIC_BC_HASH36_VALIDATED places=48 reps=%o complete=true time=%o\n",#rs,Cputime(t0);
quit;
