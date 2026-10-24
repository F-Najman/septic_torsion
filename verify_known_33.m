/* Global known C33 model and full torsion bound used in the class proof. */
R33<z33>:=PolynomialRing(Rationals());
f33:=z33^7-6*z33^6+9*z33^5-5*z33^3-z33+1;
y33:=(-3*z33^6+16*z33^5-17*z33^4-9*z33^3+11*z33^2+3*z33+4)/3;
assert IsIrreducible(f33);K33<a33>:=NumberField(f33);
function C33E(a,y)
 r:=(a^2*y-a*y+y-1)/(a^2*y-a);s:=(a*y-y+1)/(a*y);
 b:=r*s*(r-1);c:=s*(r-1);return EllipticCurve([1-c,-b,-b,0,0]);
end function;
function C33Evaluate(f,a)
 L:=Parent(a);return &+[L!Numerator(Coefficient(f,i))/(L!Denominator(Coefficient(f,i)))*a^i:i in [0..Degree(f)]];
end function;
E33:=C33E(a33,Evaluate(y33,a33));P33:=E33![0,0];
assert 33*P33 eq E33!0 and 11*P33 ne E33!0 and 3*P33 ne E33!0;
assert Degree(MinimalPolynomial(jInvariant(E33))) eq 7;
upper33:=0;
for pp33 in [5,13] do
 assert Integers()!Discriminant(f33) mod pp33 ne 0;
 targetdegree:=pp33 eq 5 select 3 else 7;
 fs33:=[fac[1]:fac in Factorization(Polynomial(GF(pp33),f33))|Degree(fac[1]) eq targetdegree];
 assert #fs33 eq 1;L33<b33>:=ext<GF(pp33)|fs33[1]>;
 EE33:=C33E(b33,C33Evaluate(y33,b33));assert Discriminant(EE33) ne 0;
 order33:=#EE33;assert order33 eq (pp33 eq 5 select 132 else 62752635);
 upper33:=GCD(upper33,order33);
 print "C33_KNOWN_GLOBAL_REDUCTION",pp33,targetdegree,order33,upper33;
end for;
assert upper33 eq 33;
print "C33_KNOWN_GLOBAL_MODEL_DONE field_degree=7 exact_point_order=33 full_torsion=33 complete=true";
