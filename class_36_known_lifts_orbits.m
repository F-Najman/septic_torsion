/* Compare the exact finite classes with the two archived global
   septimic diamond orbits. The displayed models have good reduction5. */
R36<z36>:=PolynomialRing(Rationals());
models36:=[*
 <z36^7-z36^6+2*z36^5-4*z36^4+3*z36^3+2*z36^2-2*z36+1,
  (-2*z36^6+3*z36^5-5*z36^4+11*z36^3-9*z36^2+4*z36+7)/6>,
 <z36^7-2*z36^6+4*z36^5-z36^4-8*z36^3+2*z36^2+5*z36+1,
  (-13*z36^6+36*z36^5-79*z36^4+73*z36^3+51*z36^2-65*z36-16)/2>
*];
function C36ReducePoly(f,a)
 L:=Parent(a);
 return &+[L!Numerator(Coefficient(f,i))/(L!Denominator(Coefficient(f,i)))*a^i:i in [0..Degree(f)]];
end function;
knownclasses:={G|};knownorbits:=[* *];
for idx36 in [1..2] do
 f36,y36:=Explode(models36[idx36]);assert IsIrreducible(f36);
 f5:=Polynomial(GF(5),f36);assert IsSquarefree(f5);
 D36:=DivisorGroup(C)!0;
 for fg36 in Factorization(f5) do
  g36:=fg36[1];assert fg36[2] eq 1;
  if Degree(g36) eq 1 then
   L36:=GF(5);aa36:=-Coefficient(g36,0)/Coefficient(g36,1);
  else L36<aa36>:=ext<GF(5)|g36>;end if;
  yy36:=C36ReducePoly(y36,aa36);
  pp36:=Places(C(L36)![aa36,yy36]);
  assert #pp36 eq 1 and Degree(pp36[1]) eq Degree(g36) and not IsCusp(X,pp36[1]);
  EE36:=EllipticCurve(X,pp36[1]);LL36:=LevelStructure(X,pp36[1]);
  assert Discriminant(EE36) ne 0;
  D36+:=Divisor(pp36[1]);
 end for;
 assert Degree(D36) eq 7;
 orbit36:={DiamondOperator(X,u36,D36):u36 in [1..18]|GCD(u36,36) eq 1};
 classes36:={m2(dd36-7*Divisor(c0)):dd36 in orbit36};
 print "C36_DEGREE7_KNOWN_ORBIT",idx36,[Degree(pp36):pp36 in Support(D36)],#orbit36,#classes36,Cputime(t0);
 assert #orbit36 eq 6 and #classes36 eq 6;
 assert IsEmpty(knownclasses meet classes36);
 knownclasses join:=classes36;Append(~knownorbits,orbit36);
end for;
unknown:=candidates diff knownclasses;
print "C36_DEGREE7_KNOWN_COMPARISON",#candidates,#knownclasses,"unknown",#unknown,Cputime(t0);
print "C36_DEGREE7_UNKNOWN_CLASSES",[Eltseq(gc):gc in unknown];
if IsEmpty(unknown) then
 print "C36_SEPTIMIC_COMPLETENESS_CLASSES_DONE known_orbits=2 known_classes=12 unknown=0 complete=true",Cputime(t0);
end if;
