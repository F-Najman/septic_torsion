/* Exact characteristic-zero point and full-torsion verification for
   both known septimic X1(36) points. Independent of finite class maps. */
Rglobal36<zg36>:=PolynomialRing(Rationals());
global36:=[*
 <zg36^7-zg36^6+2*zg36^5-4*zg36^4+3*zg36^3+2*zg36^2-2*zg36+1,
  (-2*zg36^6+3*zg36^5-5*zg36^4+11*zg36^3-9*zg36^2+4*zg36+7)/6>,
 <zg36^7-2*zg36^6+4*zg36^5-zg36^4-8*zg36^3+2*zg36^2+5*zg36+1,
  (-13*zg36^6+36*zg36^5-79*zg36^4+73*zg36^3+51*zg36^2-65*zg36-16)/2>
*];
function C36GlobalCurve(x,y)
 r:=(x^2*y-x*y+y-1)/(x^2*y-x);s:=(x*y-y+1)/(x*y);
 b:=r*s*(r-1);c:=s*(r-1);
 return EllipticCurve([1-c,-b,-b,0,0]);
end function;
function C36GlobalReduce(f,a)
 L:=Parent(a);
 return &+[L!Numerator(Coefficient(f,i))/(L!Denominator(Coefficient(f,i)))*a^i:i in [0..Degree(f)]];
end function;
function C36EquationValue(xcoord,ycoord)
 x:=xcoord;y:=ycoord;
 return eval Read("mdmagma/models_X1_n/FFFc36.txt");
end function;
for idxg36 in [1..2] do
 fg36,yg36:=Explode(global36[idxg36]);assert Degree(fg36) eq 7 and IsIrreducible(fg36);
 Kg36<ag36>:=NumberField(fg36);Eg36:=C36GlobalCurve(ag36,Evaluate(yg36,ag36));Pg36:=Eg36![0,0];
 assert C36EquationValue(ag36,Evaluate(yg36,ag36)) eq 0;
 assert 36*Pg36 eq Eg36!0 and 18*Pg36 ne Eg36!0 and 12*Pg36 ne Eg36!0;
 assert Degree(MinimalPolynomial(jInvariant(Eg36))) eq 7;
 boundg36:=0;
 primesg36:=idxg36 eq 1 select [29] else [5,19];
 for pg36 in primesg36 do
  assert Integers()!Discriminant(fg36) mod pg36 ne 0;
  degreeg36:=idxg36 eq 1 select 1 else (pg36 eq 5 select 3 else 2);
  gsg36:=[ff[1]:ff in Factorization(Polynomial(GF(pg36),fg36))|Degree(ff[1]) eq degreeg36];
  assert #gsg36 gt 0;gg36:=gsg36[1];
  if degreeg36 eq 1 then Lg36:=GF(pg36);aag36:=-Coefficient(gg36,0)/Coefficient(gg36,1);
  else Lg36<aag36>:=ext<GF(pg36)|gg36>;end if;
  eeg36:=C36GlobalCurve(aag36,C36GlobalReduce(yg36,aag36));
  assert Discriminant(eeg36) ne 0;
  ng36:=#eeg36;boundg36:=GCD(boundg36,ng36);
  print "C36_KNOWN_TORSION_REDUCTION",idxg36,pg36,degreeg36,ng36,boundg36;
 end for;
 assert boundg36 eq 36;
 print "C36_KNOWN_GLOBAL_MODEL_DONE",idxg36,"field_degree=7 exact_point_order=36 full_torsion=36 complete=true";
end for;
