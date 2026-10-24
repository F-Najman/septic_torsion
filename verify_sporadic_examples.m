/* Archived van Hoeij degree-seven places on X1(25), X1(33), X1(36).
   Verify irreducibility and exact point order over Q, then bound the
   complete torsion order by injective reductions at unramified odd primes. */
SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
R<x>:=PolynomialRing(Rationals());
examples:=[*
 <25,x^7+x^6-x^5-x^4+x^2+x-1,x^2+x>,
 <33,x^7-6*x^6+9*x^5-5*x^3-x+1,(-3*x^6+16*x^5-17*x^4-9*x^3+11*x^2+3*x+4)/3>,
 <36,x^7-x^6+2*x^5-4*x^4+3*x^3+2*x^2-2*x+1,(-2*x^6+3*x^5-5*x^4+11*x^3-9*x^2+4*x+7)/6>,
 <36,x^7-2*x^6+4*x^5-x^4-8*x^3+2*x^2+5*x+1,(-13*x^6+36*x^5-79*x^4+73*x^3+51*x^2-65*x-16)/2>
*];
function CurveFromXY(a,y)
 r:=(a^2*y-a*y+y-1)/(a^2*y-a);s:=(a*y-y+1)/(a*y);
 b:=r*s*(r-1);c:=s*(r-1);
 return EllipticCurve([1-c,-b,-b,0,0]);
end function;
function ReducePolynomial(f,a)
 L:=Parent(a);
 return &+[L!Numerator(Coefficient(f,i))/(L!Denominator(Coefficient(f,i)))*a^i:i in [0..Degree(f)]];
end function;
for idx in [1..#examples] do
 N,f,ypol:=Explode(examples[idx]);assert IsIrreducible(f) and Degree(f) eq 7;
 K<a>:=NumberField(f);y:=Evaluate(ypol,a);E:=CurveFromXY(a,y);P:=E![0,0];
 assert N*P eq E!0 and &and[(N div q)*P ne E!0:q in PrimeDivisors(N)];
 print "DEGREE7_SPORADIC_EXACT_POINT",idx,N,"field_degree=7","j_degree",Degree(MinimalPolynomial(jInvariant(E))),Cputime(t0);
 bound:=0;used:=0;
 for p in [3,5,7,11,13,17,19,23,29,31,37,41] do
  if Integers()!Discriminant(f) mod p eq 0 then continue;end if;
  if &or[Denominator(c) mod p eq 0:c in Coefficients(ypol)] then continue;end if;
  fac:=Factorization(Polynomial(GF(p),f));
  for fg in fac do
   g:=fg[1];dd:=Degree(g);
   if dd gt 4 and p ge 17 then continue;end if;
   if dd eq 1 then L:=GF(p);aa:=-Coefficient(g,0)/Coefficient(g,1);
   else L<aa>:=ext<GF(p)|g>;end if;
   try
    yy:=ReducePolynomial(ypol,aa);EE:=CurveFromXY(aa,yy);
    assert Discriminant(EE) ne 0;
    PP:=EE![0,0];assert N*PP eq EE!0;
    nn:=#EE;assert nn mod N eq 0;
    bound:=GCD(bound,nn);used+:=1;
    print "DEGREE7_TORSION_REDUCTION",idx,p,dd,nn,bound,Cputime(t0);
   catch er
    print "DEGREE7_BAD_DISPLAYED_MODEL",idx,p,dd;
   end try;
   if bound eq N then break;end if;
  end for;
  if bound eq N then break;end if;
 end for;
 print "DEGREE7_FULL_TORSION_BOUND",idx,N,bound,used,Cputime(t0);
 assert bound eq N;
 print "DEGREE7_SPORADIC_CERTIFIED",idx,"full_torsion",[N],"exact=true",Cputime(t0);
end for;
print "DEGREE7_SPORADIC_EXAMPLES_DONE count=4 complete=true",Cputime(t0);quit;
