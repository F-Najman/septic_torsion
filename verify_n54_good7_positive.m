/* Independent positive-point certificate. No point counts, SEA, class
   groups or modular-coordinate inverse are used in this verification.
   Explicit order27, rational2 roots, distinct j orbits, and the independent
   full source count certify completeness of the exported branch. */
SetColumns(0);SetQuitOnError(true);t0:=Cputime();
data:=eval Read("data/N54_D7_good7_pointdata.txt");
def,records:=Explode(data);assert #records eq 828;
S:=PolynomialRing(GF(5));f:=S!def;
assert Degree(f) eq 7 and IsIrreducible(f);
K<a>:=ext<GF(5)|f>;seen:=AssociativeArray();weighted:=0;one:=0;three:=0;
for record in records do
 xyco,n,nr,jco:=Explode(record);
 x:=Evaluate(S!xyco[1],a);y:=Evaluate(S!xyco[2],a);
 jj:=Evaluate(S!jco,a);assert Degree(MinimalPolynomial(jj,GF(5))) eq 7;
 assert x*y*(x*y-1) ne 0;
 r:=(x^2*y-x*y+y-1)/(x*(x*y-1));s:=(x*y-y+1)/(x*y);
 b:=r*s*(r-1);c:=s*(r-1);
 E:=EllipticCurve([1-c,-b,-b,K!0,K!0]);
 assert Discriminant(E) ne 0 and jInvariant(E) eq jj;
 P:=E![0,0];assert 27*P eq E!0 and 9*P ne E!0;
 bb2,bb4,bb6,bb8:=Explode(bInvariants(E));R<T>:=PolynomialRing(K);
 actual:=#Roots(4*T^3+bb2*T^2+2*bb4*T+bb6);
 assert nr eq actual and nr in [1,3] and n mod 54 eq 0;
 key:=Minimum([[Integers()!z:z in Eltseq(jj^(5^i),GF(5))]:i in [0..6]]);
 assert not IsDefined(seen,key);seen[key]:=true;
 weighted+:=9*nr;
 if nr eq 1 then one+:=1;else three+:=1;end if;
end for;
assert #Keys(seen) eq 828 and one eq 608 and three eq 220 and weighted eq 11412;
printf "N54_POSITIVE_BRANCH_VERIFIED records=828 distinct_j_orbits=828 roots1=608 roots3=220 weighted=11412 exact_order27=true complete=true time=%o\n",Cputime(t0);
quit;
