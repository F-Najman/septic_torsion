/* Load after the completed class_25_p2_v2.m in its retained session.
   Compare every remaining class with an actual Q-rational effective
   divisor of degree7, using the ORIGINAL diamond correspondence. */
assert #candidates eq 60 and Invariants(G) eq [227555,0];
assert Invariants(HQ) eq [227555] and IsOdd(#HQ);
RQ<zz>:=PolynomialRing(Rationals());
K3<aa>:=NumberField(zz^3-zz^2+1);assert Degree(K3) eq 3;
RK<tt>:=PolynomialRing(K3);gg:=tt^2+(aa^2-2*aa-1)*tt+aa;
assert IsIrreducible(gg);K6<bb>:=ext<K3|gg>;assert AbsoluteDegree(K6) eq 6;
xx:=K6!aa;yy:=bb;
rr:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));ss:=(xx*yy-yy+1)/(xx*yy);
bbT:=rr*ss*(rr-1);ccT:=ss*(rr-1);
E6:=EllipticCurve([1-ccT,-bbT,-bbT,K6!0,K6!0]);P6:=E6![0,0];
assert 25*P6 eq E6!0 and 5*P6 ne E6!0;
assert TateNormalForm_xy(E6,P6) eq [xx,yy];
print "C25_KNOWN_SEXTIC_GLOBAL_VERIFIED degree=6 exact_order=25";
f7Q:=zz^7+zz^6-zz^5-zz^4+zz^2+zz-1;
assert IsIrreducible(f7Q);K7<a7>:=NumberField(f7Q);y7:=a7^2+a7;
r7:=(a7^2*y7-a7*y7+y7-1)/(a7*(a7*y7-1));s7:=(a7*y7-y7+1)/(a7*y7);
b7:=r7*s7*(r7-1);c7:=s7*(r7-1);
E7:=EllipticCurve([1-c7,-b7,-b7,K7!0,K7!0]);P7:=E7![0,0];
assert 25*P7 eq E7!0 and 5*P7 ne E7!0;
assert TateNormalForm_xy(E7,P7) eq [a7,y7];
assert Degree(MinimalPolynomial(jInvariant(E7))) eq 7;
assert Integers()!Discriminant(f7Q) mod 5 ne 0;
q5:=[fac[1]:fac in Factorization(Polynomial(GF(5),f7Q))|Degree(fac[1]) eq 2];assert #q5 eq 1;
L5<a5>:=ext<GF(5)|q5[1]>;y5:=a5^2+a5;
r5:=(a5^2*y5-a5*y5+y5-1)/(a5*(a5*y5-1));s5:=(a5*y5-y5+1)/(a5*y5);
b5:=r5*s5*(r5-1);c5:=s5*(r5-1);E5:=EllipticCurve([1-c5,-b5,-b5,L5!0,L5!0]);
assert Discriminant(E5) ne 0 and #E5 eq 25;
print "C25_KNOWN_SEPTIC_GLOBAL_VERIFIED degree=7 exact_order=25 full_torsion=25 reduction_p=5 residue_degree=2 order=25";
L6:=GF(2^6);RL<z6>:=PolynomialRing(L6);
xxx:=Roots(z6^3+z6^2+1)[1][1];
rrts:=Roots(z6^2+(xxx^2+1)*z6+xxx);assert #rrts eq 2;
yyy:=rrts[1][1];pp:=Places(C(L6)![xxx,yyy]);
assert #pp eq 1 and Degree(pp[1]) eq 6 and not IsCusp(X,pp[1]);
orb6:=Seqset(DiamondOrbit(X,pp[1]));assert #orb6 eq 5 and orb6 eq Seqset(places[6]);
print "C25_KNOWN_SEXTIC_REDUCTION_ORBIT degree=6 size=5 full=true";
L7:=GF(2^7);RL7<z7>:=PolynomialRing(L7);
f7:=z7^7+z7^6+z7^5+z7^4+z7^2+z7+1;
xxx:=Roots(f7)[1][1];yyy:=xxx^2+xxx;pp:=Places(C(L7)![xxx,yyy]);
assert #pp eq 1 and Degree(pp[1]) eq 7 and not IsCusp(X,pp[1]);
orb7:=Seqset(DiamondOrbit(X,pp[1]));assert #orb7 eq 10 and orb7 eq Seqset(places[7]);
print "C25_KNOWN_SEPTIC_REDUCTION_ORBIT degree=7 size=10 full=true";
reducible:={G|};
for x in orb6 do
 for c in cs do
  if Degree(c) ne 1 then continue;end if;
  Include(~reducible,m2(Divisor(x)+Divisor(c)-7*Divisor(c0)));
 end for;
end for;
septics:={m2(Divisor(x)-7*Divisor(c0)):x in orb7};
assert #reducible eq 50 and #septics eq 10;
assert IsEmpty(reducible meet septics);
assert reducible join septics eq candidates;
print "C25_SEPTIMIC_COMPLETENESS_REDUCTION_DONE known_sextic_places=5 known_septic_places=10 reducible_classes=50 septic_classes=10 candidates=60 complete=true",Cputime(t0);
quit;
