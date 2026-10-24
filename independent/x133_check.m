SetColumns(0); load "independent/helpers7.m"; t0:=Cputime(); F2:=GF(2);
Fx<x>:=RationalFunctionField(F2); PY<Y>:=PolynomialRing(Fx);
Xpol := x^8*Y^10 - x^8*(x^9 + 2*x^8 + 4*x^7 + 6*x^6 + 9*x^5 + 12*x^4 + 16*x^3 + 20*x^2 + 25*x + 18)*Y^9
 - x*(3*x^16 + 15*x^15 + 24*x^14 + 28*x^13 + 16*x^12 - 6*x^11 - 32*x^10 - 32*x^9 + 33*x^8 + 155*x^7 + 188*x^6 + 108*x^5 + 58*x^4 + 28*x^3 + 12*x^2 + 4*x + 1)*Y^8
 - x*(6*x^16 + 48*x^15 + 144*x^14 + 226*x^13 + 225*x^12 + 39*x^11 - 402*x^10 - 1116*x^9 - 1854*x^8 - 2071*x^7 - 1536*x^6 - 852*x^5 - 476*x^4 - 249*x^3 - 114*x^2 - 40*x - 9)*Y^7
 - 2*x*(x + 1)*(5*x^15 + 50*x^14 + 210*x^13 + 522*x^12 + 981*x^11 + 1599*x^10 + 2231*x^9 + 2597*x^8 + 2292*x^7 + 1476*x^6 + 774*x^5 + 516*x^4 + 332*x^3 + 182*x^2 + 72*x + 18)*Y^6
 + x*(x + 1)*(x^18 + 9*x^17 + 51*x^16 + 204*x^15 + 591*x^14 + 1281*x^13 + 2295*x^12 + 3747*x^11 + 5451*x^10 + 6421*x^9 + 5607*x^8 + 3279*x^7 + 1568*x^6 + 1620*x^5 + 2142*x^4 + 1775*x^3 + 1026*x^2 + 396*x + 84)*Y^5
 + 2*x*(x + 1)*(5*x^15 + 25*x^14 + 35*x^13 - 67*x^12 - 278*x^11 - 235*x^10 + 431*x^9 + 1297*x^8 + 1193*x^7 - 405*x^6 - 2241*x^5 - 2775*x^4 - 2040*x^3 - 1047*x^2 - 357*x - 63)*Y^4
 + (x + 1)*(6*x^16 + 48*x^15 + 144*x^14 + 110*x^13 - 529*x^12 - 1863*x^11 - 2569*x^10 - 451*x^9 + 4950*x^8 + 10657*x^7 + 12679*x^6 + 10287*x^5 + 6104*x^4 + 2675*x^3 + 783*x^2 + 113*x - 1)*Y^3
 + (x + 1)*(3*x^16 + 33*x^15 + 159*x^14 + 413*x^13 + 471*x^12 - 621*x^11 - 3915*x^10 - 9165*x^9 - 14088*x^8 - 15812*x^7 - 13476*x^6 - 8904*x^5 - 4576*x^4 - 1752*x^3 - 432*x^2 - 44*x + 3)*Y^2
 + (x + 1)^8*(x^9 + 7*x^8 + 24*x^7 + 50*x^6 + 71*x^5 + 69*x^4 + 50*x^3 + 20*x^2 + 15*x - 3)*Y + (x + 1)^8;
assert IsIrreducible(Xpol); F<y>:=FunctionField(Xpol);
q:=(1+x)/x; t:=1/((1-y)*(1+x));
E:=EllipticCurve([(1+q)*t+(2-q),0,q*t^2+(1-q)*t,0,0]); P0:=E![-t,t^2];
j:=jInvariant(E); printf "genus %o, deg j %o\n", Genus(F), Degree(j); assert Genus(F) eq 21 and Degree(j) eq 480;
Q3:=3*P0; Q11:=11*P0; assert Q3 ne E!0 and Q11 ne E!0 and 3*Q11 eq E!0; print "base point has exact order 33";
b,c:=tate(E,P0);
cusps:=[Q: Q in Support(Divisor(j)) | Valuation(j,Q) lt 0];
printf "cusp degrees %o\n", Sort([Degree(Q): Q in cusps]);
nonc:=[[Q: Q in Places(F,d) | Valuation(j,Q) ge 0]: d in [1..7]];
printf "noncuspidal place counts %o (%o s)\n", [#L: L in nonc], Cputime(t0);
Cl,m:=ClassGroup(F); printf "class group %o (%o s)\n", Invariants(Cl), Cputime(t0);
c0:=[Q: Q in cusps | Degree(Q) eq 1][1]; cl:=func<D | (D - Degree(D)*c0) @@ m>;
five:=[Q: Q in cusps | Degree(Q) eq 5]; assert #five eq 2;
H:=sub<Cl | [cl(1*Q): Q in cusps | Degree(Q) ne 5] cat [cl(five[1]+five[2])]>;
printf "H (reductions of Q-rational cuspidal divisors): %o, order %o, DEHMZ order %o\n", Invariants(H), #H, 5*42373650;
L7:=[]; effdivs(cusps,1,7,DivisorGroup(F)!0,~L7);
liftable:=func<D | Valuation(D,five[1]) eq 0 and Valuation(D,five[2]) eq 0>;
qcl:={cl(D): D in L7 | liftable(D)};
badcusp:=[D: D in L7 | not liftable(D) and cl(D) in H and cl(D) notin qcl];
printf "effective cuspidal deg 7: %o (liftable %o, classes %o); nonliftable with class in H outside rational cuspidal classes: %o\n", #L7, #[D: D in L7 | liftable(D)], #qcl, #badcusp;
L1:=[]; effdivs(cusps,1,1,DivisorGroup(F)!0,~L1); L2:=[]; effdivs(cusps,1,2,DivisorGroup(F)!0,~L2);
in5:=#[1: Q in nonc[5], C in L2 | cl(1*Q+C) in H]; in6:=#[1: Q in nonc[6], C in L1 | cl(1*Q+C) in H];
in7:=[Q: Q in nonc[7] | cl(1*Q) in H and cl(1*Q) notin qcl];
printf "P5+C2 in H: %o of %o; P6+C1 in H: %o of %o; P7 in H (not rational cuspidal class): %o of %o\n", in5, #nonc[5]*#L2, in6, #nonc[6]*#L1, #in7, #nonc[7];
Qz<Z>:=PolynomialRing(Rationals()); K7<a>:=NumberField(Z^7-6*Z^6+9*Z^5-5*Z^3-Z+1);
yy:=(-3*a^6+16*a^5-17*a^4-9*a^3+11*a^2+3*a+4)/3;
kr:=known_reductions(K7, a, yy, 33, 2);
printf "residue degrees at 2 of the known point: %o\n", [ [t[1]: t in L][1] : L in kr];
cov:={}; for L in kr do for tr in L do if tr[1] eq 7 then h:=matchplace(nonc[7],tr,b,c); assert #h eq 1; Include(~cov,h[1]); end if; end for; end for;
printf "septic places covered by translates of the known point: %o; all P7 classes in H covered: %o\n", #cov, Set(in7) subset cov;
printf "X133_CHECK_DONE (%o s)\n", Cputime(t0); quit;
