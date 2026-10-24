SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_diamond_hash_m2.m");
Attach("prime_degree_place.m");
X:=MDX11(2,28,GF(3));cs:=Cusps(X);ps:=[NoncuspidalPlaces(X,d):d in [1..7]];
assert [#s:s in ps] eq [0,0,0,0,48,64,312];
rs:=PlacesUpToDiamondM2HashDirect(X,ps[7]);assert #rs eq 52;
K:=Parent(Evaluate(X`_coordinates[1],rs[1]));assert #K eq 3^7;
records:=[* *];
for i in [1..52] do
 bc:=MDM2CanonicalPair(X,rs[i],K);
 pp:=MDPlaceCommonZerosPrimeDegree(X`_coordinates,bc);
 assert pp eq rs[i] and Degree(pp) eq 7;
 Append(~records,[[Integers()!c:c in Eltseq(z,GF(3))]:z in bc]);
 print "N28_D7_FIXED_INPUT",i,"roundtrip=true",Cputime(t0);
end for;
def:=[Integers()!c:c in Coefficients(DefiningPolynomial(K,GF(3)))];
Write("data/N28_D7_hecke_input.txt",Sprint([*def,records*]) : Overwrite:=true);
print "N28_D7_FIXED_INPUT_DONE records=52 complete=true",Cputime(t0);quit;
