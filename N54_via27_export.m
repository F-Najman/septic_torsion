/* Degree7 necessary two-prime class test on X1(27).
   The separate semistable-model argument forces some good reduction at2.
   At5 an irreducible 2-division cubic permits only multiplicity divisible3.
   Cuspidal complements are all included. */
SetColumns(0);SetMemoryLimit(16*10^9);SetQuitOnError(true);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("prime_degree_place.m");
load "load_n54_good7.m";
function CuspDivisors(cs,d,C)
 if d eq 0 then return [DivisorGroup(C)!0];end if;
 if #cs eq 0 then return [];end if;
 P:=cs[1];rest:=cs[2..#cs];out:=[];
 for k in [0..d div Degree(P)] do
  out cat:=[k*Divisor(P)+a:a in CuspDivisors(rest,d-k*Degree(P),C)];
 end for;
 return out;
end function;

X2:=MDX1(27,GF(2));C2:=Curve(X2);
os2:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in CuspOrbitsQ(X2)];
sig2:=[CuspSignature(X2,Support(o)[1]):o in os2];
assert #Seqset(sig2) eq #os2;
c02:=Support([o:o in os2|Degree(o) eq 1][1])[1];
baseSig:=CuspSignature(X2,c02);
G2,a2,b2:=ClassGroupRationalReduction(C2,c02:Proof:=true);
qc2:=[b2(o-Degree(o)*Divisor(c02)):o in os2];
H2:=sub<G2|qc2>;assert Invariants(H2) eq [3,3,52497];
print "N54_D7_F2_GROUP",Invariants(G2),Invariants(H2),Cputime(t0);
ps2:=[NoncuspidalPlaces(X2,d):d in [1..7]];
assert [#p:p in ps2[1..6]] eq [0,0,0,0,0,9];
candidate2:={G2|};
for d in [6,7] do
 for p in ps2[d] do
  for cc in CuspDivisors(Cusps(X2),7-d,C2) do
   a:=b2(Divisor(p)+cc-7*Divisor(c02));
   if a in H2 then Include(~candidate2,a);end if;
  end for;
 end for;
end for;
print "N54_D7_F2_CANDIDATES",[#p:p in ps2],#candidate2,Cputime(t0);

X5:=MDX1(27,GF(5));C5:=Curve(X5);
os5raw:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in CuspOrbitsQ(X5)];
sig5:=[CuspSignature(X5,Support(o)[1]):o in os5raw];
assert Seqset(sig5) eq Seqset(sig2) and #Seqset(sig5) eq #os5raw;
os5:=[os5raw[Index(sig5,s)]:s in sig2];
c05:=Support(os5[Index(sig2,baseSig)])[1];assert Degree(c05) eq 1;
G5,a5,b5:=ClassGroupRationalReduction(C5,c05:Proof:=true);
qc5:=[b5(o-Degree(o)*Divisor(c05)):o in os5];
H5:=sub<G5|qc5>;assert Invariants(H5) eq [3,3,52497];
L:=FreeAbelianGroup(#os2);
h2:=hom<L->G2|qc2>;h5:=hom<L->G5|qc5>;
assert Image(h2) eq H2 and Image(h5) eq H5;
assert Kernel(h2) eq Kernel(h5);
print "N54_D7_GLOBAL_CUSP_RELATIONS_MATCH",#os2;
print "N54_D7_F5_GROUP",Invariants(G5),Invariants(H5),Cputime(t0);
out54:=Open("data/N54_class_arithmetic_data.py","w");
fprintf out54,"relations2 = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G2))];
fprintf out54,"relations5 = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G5))];
fprintf out54,"cusp2 = %o\n",[Eltseq(a):a in qc2];
fprintf out54,"cusp5 = %o\n",[Eltseq(a):a in qc5];
fprintf out54,"candidate2 = %o\n",[Eltseq(a):a in candidate2];
fprintf out54,"tests = [\n";

cs5:=Cusps(X5);comps:=AssociativeArray();
for k in [0..4] do
 comps[k]:={b5(c-k*Divisor(c05)):c in CuspDivisors(cs5,k,C5)};
 print "N54_D7_CUSP_CLASSES",k,#comps[k],Cputime(t0);
end for;
allps:=[[]:d in [1..7]];
for d in [1..6] do
 print "N54_FAST_LOWER_PLACE_START",d,Cputime(t0);
 allps[d]:=NoncuspidalPlaces(X5,d);
 print "N54_FAST_LOWER_PLACE_DONE",d,#allps[d],Cputime(t0);
end for;
assert #allps[1] eq 0;
good:=[[]:d in [1..7]];inert2:=[];rootsum:=[0:d in [1..7]];rootone3:=0;
for d in [2..6] do
 print "N54_FAST_2ROOT_START",d,Cputime(t0);
 for p in allps[d] do
  E:=EllipticCurve(X5,p);
  bb2,bb4,bb6,bb8:=Explode(bInvariants(E));R2<x2>:=PolynomialRing(BaseRing(E));
  f2:=4*x2^3+bb2*x2^2+2*bb4*x2+bb6;
  nr:=#Roots(f2);assert nr in [0,1,3];rootsum[d]+:=nr;
  if d eq 3 and nr eq 1 then rootone3+:=1;end if;
  has2:=nr gt 0;
  if has2 then Append(~good[d],p);
  elif d eq 2 then Append(~inert2,p);end if;
 end for;
 print "N54_FAST_2ROOT_DONE",d,#good[d],Cputime(t0);
end for;
assert #good[2] eq 0;
assert rootsum[3..5] eq [54,261,612];
assert rootsum[6]+#inert2+rootone3 eq 3033;
print "N54_D7_LOWER_SOURCE_COUNTS",rootsum[3..6],#inert2,rootone3;
good[7]:=N54LoadGood7(X5,"data/N54_D7_good7_pointdata.txt");
print "N54_D7_F5_PLACE_COUNTS",[#p:p in allps],[#p:p in good],#inert2,Cputime(t0);
checked:=0;members:=0;passing:=0;allhits:=[];
procedure TestClass(a,e,desc,~checked,~members,~passing,~allhits)
 for cc in comps[7-e] do
  checked+:=1;q:=a+cc;is_member:=0;is_pass:=0;im2:=[];
  if q in H5 then
   members+:=1;is_member:=1;v:=q@@h5;assert h5(v) eq q;im2:=Eltseq(h2(v));
   if h2(v) in candidate2 then
    passing+:=1;is_pass:=1;Append(~allhits,<desc,Eltseq(v)>);
    print "N54_D7_TWO_PRIME_HIT",desc,Eltseq(v);
   end if;
  end if;
  fprintf out54,"[%o,%o,%o,%o,%o],\n",e,Eltseq(q),is_member,is_pass,im2;
 end for;
end procedure;
for d in [3..7] do
 rs:=d eq 7 select good[7] else PlacesUpToDiamondFast(X5,good[d]);
 print "N54_D7_SINGLE_START",d,#good[d],#rs,Cputime(t0);
 for i in [1..#rs] do
  a:=b5(Divisor(rs[i])-d*Divisor(c05));
  TestClass(a,d,<d,i>,~checked,~members,~passing,~allhits);
  if i mod 10 eq 0 then print "N54_D7_SINGLE_PROGRESS",d,i,#rs,checked,members,passing,Cputime(t0);end if;
 end for;
 print "N54_D7_SINGLE_DONE",d,#rs,checked,members,passing,Cputime(t0);
end for;
for p in PlacesUpToDiamondFast(X5,inert2) do
 a:=b5(3*Divisor(p)-6*Divisor(c05));
 TestClass(a,6,<2,3>,~checked,~members,~passing,~allhits);
end for;
/* In mixed support reduce only the first place by diamonds. */
q3:=[b5(Divisor(p)-3*Divisor(c05)):p in good[3]];
q4:=[b5(Divisor(p)-4*Divisor(c05)):p in good[4]];
for p in PlacesUpToDiamondFast(X5,good[3]) do
 a:=b5(Divisor(p)-3*Divisor(c05));
 for b in q3 do TestClass(a+b,6,<3,3>,~checked,~members,~passing,~allhits);end for;
 for b in q4 do TestClass(a+b,7,<3,4>,~checked,~members,~passing,~allhits);end for;
end for;
printf "N54_D7_TWO_PRIME_DONE checked=%o members=%o passing=%o complete=true time=%o\n",checked,members,passing,Cputime(t0);
fprintf out54,"]\ncomplete = True\n";delete out54;
Write("data/N54_two_prime_export_hits.txt",Sprint(allhits) : Overwrite:=true);
print allhits;
quit;
