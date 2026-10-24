SetColumns(0);
Qx<x>:=PolynomialRing(Rationals());
data:=[ <25, x^7+x^6-x^5-x^4+x^2+x-1, func<a| a^2+a>>,
        <33, x^7-6*x^6+9*x^5-5*x^3-x+1, func<a| (-3*a^6+16*a^5-17*a^4-9*a^3+11*a^2+3*a+4)/3>>,
        <36, x^7-x^6+2*x^5-4*x^4+3*x^3+2*x^2-2*x+1, func<a| (-2*a^6+3*a^5-5*a^4+11*a^3-9*a^2+4*a+7)/6>>,
        <36, x^7-2*x^6+4*x^5-x^4-8*x^3+2*x^2+5*x+1, func<a| (-13*a^6+36*a^5-79*a^4+73*a^3+51*a^2-65*a-16)/2>> ];
for d in data do
  N:=d[1]; f:=d[2]; assert IsIrreducible(f);
  K<a>:=NumberField(f); y:=d[3](a);
  r:=(a^2*y-a*y+y-1)/(a^2*y-a); s:=(a*y-y+1)/(a*y);
  b:=r*s*(r-1); c:=s*(r-1);
  E:=EllipticCurve([1-c,-b,-b,0,0]); P:=E![0,0];
  assert Order(P) eq N;
  OK:=MaximalOrder(K);
  T,mT:=TorsionSubgroup(E);
  printf "N=%o disc(K)=%o j-degree=%o torsion=%o\n", N, Discriminant(OK), Degree(MinimalPolynomial(jInvariant(E))), Invariants(T);
end for;
print "EXAMPLES7_DONE";
quit;
