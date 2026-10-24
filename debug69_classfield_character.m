SetColumns(0);SetQuitOnError(false);SetMemoryLimit(16*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
if assigned N then N:=StringToInteger(N);else N:=57;end if;
assert N in [57,69];m:=N eq 55 select 5 else 3;ell:=N div m;h:=(ell-1) div 2;tstart:=Cputime();
CQ:=SmallModularCurve(ell);jraw:=jInvariant(CQ,ell);wraw:=DefiningPolynomials(AtkinLehnerInvolution(CQ,ell,ell));
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq (ell eq 23 select 2 else 1);
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);ww:=[RC!e:e in wraw];
jell:=FunctionField(CC)!(jn/jd);jNell:=FunctionField(CC)!(Evaluate(jn,ww)/Evaluate(jd,ww));
K,mK:=AlgorithmicFunctionField(FunctionField(CC));jK:=mK(jell);
jdegree:=-&+[Valuation(jK,P)*Degree(P):P in Poles(jK)];assert jdegree eq ell+1;
printf "DATABASE_ELL_MODEL N=%o ell=%o genus=%o good_reduction=true j_degree=%o time=%o\n",N,ell,Genus(K),jdegree,Cputime(tstart);
R<v>:=PolynomialRing(K);f:=N eq 55 select v^6+jK*v+jK else v^4+v^3+1/jK;assert IsIrreducible(f);
AF<t3>:=FunctionField(f);gg:=Genus(AF);assert gg eq (N eq 69 select 7 else 5);
O:=MaximalOrderFinite(AF);OK:=MaximalOrderFinite(K);
numY,denY:=IntegralSplit(K.1,OK);numJ,denJ:=IntegralSplit(jK,OK);
IX:=AF!BaseField(K).1;IY:=AF!K!numY;multT:=N eq 55 select K!denJ else K!numJ;IT:=(AF!multT)*t3;
integral_coordinates:=[O!IX,O!IY,O!IT];
fx:=(K!BaseField(K).1)@@mK;fy:=(K!numY)@@mK;fJ:=multT@@mK;
printf "X0_3ELL_MODEL N=%o genus=%o degree_over_X0ell=4 integral_coordinates_ready=true time=%o\n",N,gg,Cputime(tstart);
function ThreeParameter(E,R3)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x:=R3[1];y:=R3[2];
 aa3:=2*y+a3+a1*x;aa2:=3*x+a2;aa4:=3*x^2+2*x*a2+a4-a1*y;
 s:=aa4/aa3;aone:=a1+2*s;assert aa2-s^2-a1*s eq 0;
 tt:=aa3/aone^3;assert 1/(tt^3*(tt+1)) eq jInvariant(E);return tt;
end function;
function SmallParameter(E,R)
 if N eq 55 then
  b,c:=Explode(TateNormalForm_bc(E,R));assert b eq c;u:=b+1/b;assert u^6/(u+1) eq jInvariant(E);return u;
 else return ThreeParameter(E,R);end if;
end function;
function EllPoint(E,Rl)
 L:=BaseRing(E);PL<X>:=PolynomialRing(L);polyker:=&*[X-(k*Rl)[1]:k in [1..(ell-1) div 2]];
 Et,phi:=IsogenyFromKernel(E,polyker);jj:=jInvariant(E);jjN:=jInvariant(Et);
 A:=Evaluate(jn,[X,PL!0,PL!1]);BB:=Evaluate(jn,[X,PL!1,PL!1])-A;DD:=Evaluate(jd,[X,PL!0,PL!1]);
 assert Evaluate(jd,[X,PL!1,PL!1]) eq DD;
 ce:=DefiningEquation(CC);c0:=Evaluate(ce,[X,PL!0,PL!1]);c2:=PL!1;
 if ell eq 19 then c1:=PL!1;else c1:=X^3+X+1;end if;
 assert Evaluate(ce,[X,PL!1,PL!1]) eq c0+c1+c2;
 YY:=jj*DD-A;res:=c2*YY^2+c1*YY*BB+c0*BB^2;
 out:=[];
 for rt in Roots(res) do
  xx:=rt[1];bb:=Evaluate(BB,xx);
  if bb ne 0 then ys:=[Evaluate(YY,xx)/bb];else ys:=[zz[1]:zz in Roots(c2*X^2+Evaluate(c1,xx)*X+Evaluate(c0,xx))];end if;
  for yy in ys do
   pt:=CC(L)![xx,yy,L!1];
   if Evaluate(jd,Eltseq(pt)) ne 0 and Evaluate(jell,pt) eq jj and Evaluate(jNell,pt) eq jjN then Append(~out,pt);end if;
  end for;
 end for;
 assert #out eq 1;return out[1];
end function;
function Coordinates0(E,R3,Rl)
 pt:=EllPoint(E,Rl);return <Evaluate(fx,pt),Evaluate(fy,pt),Evaluate(fJ,pt)*SmallParameter(E,R3)>;
end function;
F<a>:=GF(128);lookup:=AssociativeArray();models:=AssociativeArray();
function Key0(r)
 return Minimum([&cat[Eltseq(s^(2^i),GF(2)):s in r]:i in [0..6]]);
end function;
count:=0;
for j in F do
 E0:=EllipticCurveFromjInvariant(j);
 for E in Twists(E0) do
  if #E mod N ne 0 then continue;end if;assert #E eq 2*N;count+:=1;
  A,mp:=AbelianGroup(E);assert Invariants(A) eq [2*N];PN:=2*mp(A.1);
  R3:=ell*PN;Rl:=m*PN;r:=Coordinates0(E,R3,Rl);key:=Key0(r);
  lookup[key]:=r;models[key]:=<E,R3,Rl>;
 end for;
end for;
assert count eq (N eq 55 select 7 else (N eq 57 select 14 else 21)) and #Keys(lookup) eq (N eq 55 select 1 else (N eq 57 select 2 else 3));
printf "X0_3ELL_NC7 N=%o curves=%o closed_image_points=%o time=%o\n",N,count,#Keys(lookup),Cputime(tstart);
function CoordPlace(r)
 L:=Parent(r[1]);degrees:=[Degree(MinimalPolynomial(s,GF(2))):s in r];d:=LCM(degrees);
 candidates:=[<IX,r[1]>,<IY,r[2]>,<IT,r[3]>,<IX+IY,r[1]+r[2]>,<IX+IT,r[1]+r[3]>,<IY+IT,r[2]+r[3]>,<IX+IY+IT,r[1]+r[2]+r[3]>];
 found:=false;
 for q in candidates do
  mp:=MinimalPolynomial(q[2],GF(2));if Degree(mp) eq d then gen:=q[1];val:=q[2];found:=true;break;end if;
 end for;
 assert found;
 MB:=Matrix(GF(2),[Eltseq(val^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 gens:=[O!Evaluate(mp,gen)];
 for i in [1..3] do
  xv:=Vector(GF(2),Eltseq(r[i],GF(2)));assert Rank(VerticalJoin(MB,Matrix(xv))) eq d;
  co:=Solution(MB,xv);poly:=&+[AF!co[k+1]*gen^k:k in [0..d-1]];
  Append(~gens,O!([IX,IY,IT][i]-poly));
 end for;
 II:=ideal<O|gens>;assert IsPrime(II);P:=Place(II);assert Degree(P) eq d;return P;
end function;

cs:=Poles(AF!jK);assert #cs eq 4 and &and[Degree(c) eq 1:c in cs];
Dzero:=DivisorGroup(AF)!0;RG,rm:=RayClassGroup(Dzero);Q,q:=quo<RG|[2*g:g in Generators(RG)]>;
assert Invariants(RG) eq [22,44,0] and Invariants(Q) eq [2,2,2];
cols:=[Vector(GF(2),Eltseq(q(Divisor(c)@@rm))):c in cs];
expected:=[GF(2)|(-Valuation(AF!jK,c) in [1,23]) select 1 else 0:c in cs];
chars:=[v:v in VectorSpace(GF(2),3)|&and[InnerProduct(v,cols[i]) eq expected[i]:i in [1..4]]];
print "C69_COVER_CUSP_CHARACTERS",#chars,chars;assert #chars eq 2;
ps:=[CoordPlace(lookup[key]):key in Sort(Setseq(Keys(lookup)))];
pcs:=[Vector(GF(2),Eltseq(q(Divisor(P)@@rm))):P in ps];
chars:=[v:v in chars|&and[InnerProduct(v,cl) eq 0:cl in pcs]];
assert #chars eq 1;chi:=chars[1];
print "C69_COVER_IDENTIFIED",chi,"cusp_widths",[-Valuation(AF!jK,c):c in cs],"cusp_values",expected,"nc7_classes",pcs,Cputime(tstart);

C2:=AbelianGroup([2]);hh:=hom<RG->C2|[C2![Integers()!InnerProduct(chi,Vector(GF(2),Eltseq(q(RG.i))))]:i in [1..Ngens(RG)]]>;
assert &and[hh(Divisor(cs[i])@@rm) eq C2![Integers()!expected[i]]:i in [1..#cs]];
assert &and[hh(Divisor(P)@@rm) eq C2!0:P in ps];
AA:=AbelianExtension(Dzero,Kernel(hh));CYns:=FunctionField(AA);pp:=DefiningPolynomials(CYns);assert #pp eq 1 and Degree(pp[1]) eq 2;CY:=FunctionField(pp[1]);assert Genus(CY) eq 13;
print "C69_FRESH_CLASS_FIELD_READY",Genus(CY),Cputime(tstart);

print "DEBUG_POLY_TYPE",Type(pp[1]),CoefficientRing(Parent(pp[1])) eq AF;
print "DEBUG_BASE_FIELDS",BaseField(CY) eq AF,BaseField(CYns) eq AF;
print "DEBUG_ROOTS",Evaluate(pp[1],CY.1),Evaluate(pp[1],CYns.1);
for c in cs do
 width:=-Valuation(AF!jK,c);pred:=hh(Divisor(c)@@rm);
 print "DEBUG_CUSP_START",width,"expected",pred,"abstract_type",DecompositionType(AA,c),Cputime(tstart);
 dd:=Decomposition(CY,c);print "DEBUG_CUSP_ACTUAL",width,#dd,[Degree(a):a in dd],Cputime(tstart);
end for;
for i in [1..Ngens(RG)] do
 Dg:=rm(RG.i);fr:=0;details:=[];
 for B in Support(Dg) do
  dd:=Decomposition(CY,B);Append(~details,<Degree(B),Valuation(Dg,B),#dd,[Degree(a):a in dd]>);
  if #dd eq 1 then fr+:=Valuation(Dg,B);end if;
 end for;
 print "DEBUG_GENERATOR",i,"expected",hh(RG.i),"observed",fr mod 2,"details",details,Cputime(tstart);
end for;
print "C69_CLASS_FIELD_DEBUG_DONE",Cputime(tstart);
