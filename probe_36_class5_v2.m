/* Exact bounded first probe for septimic points on X1(36).
   This source is independent of the degree-six computations. */
SetColumns(0);SetMemoryLimit(16*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
fp:=FrobeniusPolynomial(JOne(36),5);jf5:=Integers()!Evaluate(fp,1);
assert jf5 eq 5061852020736;
print "C36_F5_JORDER",jf5,Factorization(jf5),Cputime(t0);
X:=MDX1(36,GF(5));C:=Curve(X);cs:=Cusps(X);os:=CuspOrbitsQ(X);
assert Genus(X) eq 17;
assert Sort([Degree(o):o in os]) eq [1:i in [1..9]] cat [2:i in [1..6]] cat [3,4,6,6];
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
print "C36_F5_SETUP",[Degree(c):c in cs],[Degree(o):o in os],Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[n:n in Invariants(G)|n ne 0] eq jf5;
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
qc:=[m2(o-Degree(o)*Divisor(c0)):o in os];
Hall:=sub<G|cc>;HQ:=sub<G|qc>;
assert Invariants(HQ) eq [12,252,7812];
print "C36_F5_CLASS_GROUPS",Invariants(G),Invariants(Hall),Invariants(HQ),Cputime(t0);
places:=[NoncuspidalPlaces(X,d):d in [1..7]];
assert [#ss:ss in places][1..6] eq [0,0,48,144,660,2856];
print "C36_F5_NONCUSP",[#ss:ss in places],Cputime(t0);
print "C36_F5_PROBE_DONE complete=true",Cputime(t0);
