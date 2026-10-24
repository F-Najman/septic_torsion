/* Exact necessary-condition sieve: all effective degree7 divisors in the rational
   cuspidal source-image bound on the genus9 A3 quotient. */
SetColumns(0);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("class_rational_reduction.m");
Attach("fast_diamond_hash_m2.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("fast_hecke_cyclic_commonzero_v1.m");
load "a3_m2_30_model.m";
G,_,phi:=ClassGroup(A3K:ReductionDivisor:=Divisor(A3c0),Proof:=true);
hgens:=[phi(OO-Degree(OO)*Divisor(A3c0)):OO in A3orbits];H:=sub<G|hgens>;
assert Invariants(H) eq [2,2,2,2,2,8] and #H eq 256;
print "A3_D7_RR_GROUP",Invariants(G),Invariants(H),Cputime(t0);
Free:=FreeAbelianGroup(#A3orbits);hm:=hom<Free->G|hgens>;assert Image(hm) eq H;
KB:=Matrix(Integers(),[Eltseq(Free!a):a in Generators(Kernel(hm))]);
Lat:=LatticeWithBasis(LLL(KB));
classes:=[OO-Degree(OO)*Divisor(A3c0):OO in A3orbits];
hist:=AssociativeArray();records:=[* *];targetdivisors:=[];noncusp:=0;tested:=0;
for hh in H do
 tested+:=1;vv:=Eltseq((G!hh)@@hm);ww:=ClosestVector(Lat,Vector(Rationals(),vv));
 vv:=[vv[i]-Integers()!ww[i]:i in [1..#vv]];assert hm(Free!vv) eq G!hh;
 DD:=7*Divisor(A3c0)+&+[vv[i]*classes[i]:i in [1..#vv]];assert Degree(DD) eq 7;
 V,vm:=RiemannRochSpace(DD);dim:=Dimension(V);
 if not IsDefined(hist,dim) then hist[dim]:=0;end if;hist[dim]+:=1;
 patterns:=[* *];
 if dim gt 0 then
  for pt in Points(ProjectiveSpace(GF(7),dim-1)) do
   ff:=vm(V!Eltseq(pt));TT:=DD+Divisor(ff);assert IsEffective(TT) and Degree(TT) eq 7;
   ss,mm:=Support(TT);nc:=&or[q notin A3cs:q in ss];
   Append(~patterns,[*[Degree(q):q in ss],mm,nc*]);Append(~targetdivisors,TT);
   if nc then noncusp+:=1;end if;
  end for;
  print "A3_D7_RR_CLASS",tested,dim,#patterns,Cputime(t0);
 end if;
 Append(~records,[*tested,Eltseq(G!hh),dim,patterns*]);
end for;
assert tested eq 256 and #Seqset(targetdivisors) eq #targetdivisors;
print "A3_D7_RR_HISTOGRAM",[<k,hist[k]>:k in Sort(Setseq(Keys(hist)))];
Write("data/A3_D7_RR_records_v2.txt",Sprint(records) : Overwrite:=true);
printf "A3_D7_RR_DONE classes=%o divisors=%o with_noncuspidal_support=%o complete=true time=%o\n",tested,#targetdivisors,noncusp,Cputime(t0);
Y:=MDX1(30,GF(7));CY:=Curve(Y);assert Genus(Y) eq 9;
osY:=CuspOrbitsQ(Y);c0Y:=Support([o:o in osY|Degree(o) eq 1][1])[1];
GY,m1Y,m2Y:=ClassGroupRationalReduction(CY,c0Y:Proof:=true);
assert Invariants(GY) eq [4,4,4,136,16320,0];
hgensY:=[m2Y(o-Degree(o)*Divisor(c0Y)):o in osY];
HY:=sub<GY|hgensY>;assert Invariants(HY) eq [4,8160];QY,qY:=quo<GY|HY>;
csY:=Cusps(Y);dsY:=[Degree(x):x in csY];
ccY:=[m2Y(Divisor(x)-Degree(x)*Divisor(c0Y)):x in csY];
procedure AddCuspY(remaining,start,current,~out)
 if remaining eq 0 then Include(~out,current);return;end if;
 for i in [start..#dsY] do
  if dsY[i] le remaining then AddCuspY(remaining-dsY[i],i,current+qY(ccY[i]),~out);end if;
 end for;
end procedure;
complementsY:=AssociativeArray();
for k in [0..5] do
 out:={QY!0};Exclude(~out,QY!0);AddCuspY(k,1,QY!0,~out);complementsY[k]:={-a:a in out};
 if k le 4 then assert #out eq [1,9,41,67,68][k+1];end if;
 print "M230_D7_Y_COMPLEMENT",k,#out;
end for;
print "M230_D7_Y_GROUP",Invariants(GY),Invariants(HY),Cputime(t0);
// One entry per actual A3 closed place: the E15 fibre and its cubic2division polynomial.
infoCache:=AssociativeArray();mapCache:=AssociativeArray();localRecords:=[* *];
function FibreInfo(q)
 dd:=Degree(q);assert dd in [2..7] and q notin A3cs;
 ai:=[Evaluate(A3K!a,q):a in aInvariants(A3E15)];
 pv:=[Evaluate(A3K!a,q):a in Eltseq(A3P15)];
 vals:=MDM2ToCommonField(ai cat pv cat [Evaluate(A3root,q)],A3fields[dd]);
 ai:=vals[1..5];pv:=vals[6..8];root:=vals[9];K:=A3fields[dd];
 assert root^2 eq Discriminant(EllipticCurve(ai));
 a1,a2,a3,a4,a6:=Explode(ai);
 // Correct square completion. No generic a1!=0 M2 inverse is used.
 shortai:=[K!0,a2+a1^2/4,K!0,a4+a1*a3/2,a6+a3^2/4];
 EE:=EllipticCurve(shortai);PP:=EE![pv[1],pv[2]+(a1*pv[1]+a3)/2];
 assert 15*PP eq EE!0 and 5*PP ne EE!0 and 3*PP ne EE!0;
 assert Discriminant(EE) eq root^2 and root ne 0;
 R<T>:=PolynomialRing(K);ff:=T^3+shortai[2]*T^2+shortai[4]*T+shortai[5];
 fac:=Factorization(ff);ds:=Sort([Degree(a[1]):a in fac]);assert ds in [[1,1,1],[3]];
 return [*dd,EE,PP,fac,ds eq [1,1,1]*];
end function;
function CyclicClass(EE,PP,dd)
 assert 30*PP eq EE!0 and 15*PP ne EE!0 and 10*PP ne EE!0 and 6*PP ne EE!0;
 yy:=ModuliPointCyclicCommonZeros(Y,EE,rec<recformat<P:PtEll>|P:=PP>);
 dy:=Degree(yy);assert dd mod dy eq 0 and not IsCusp(Y,yy);
 return (dd div dy)*m2Y(Divisor(yy)-dy*Divisor(c0Y)),dy;
end function;
procedure LocalOptions(q,mul,~infoCache,~mapCache,~localRecords,~options)
 if not IsDefined(infoCache,q) then infoCache[q]:=FibreInfo(q);end if;
 dd,EE,PP,fac,split:=Explode(infoCache[q]);
 options:={GY!0};Exclude(~options,GY!0);
 if not split and mul mod 3 ne 0 then return;end if;
 if not IsDefined(mapCache,q) then
  cls:=[];imageDegrees:=[];
  if split then
   for f in fac do
    xx:=-Coefficient(f[1],0)/Coefficient(f[1],1);TT:=EE![xx,0];QQ:=8*PP+TT;
    assert 2*QQ eq PP;
    cl,dy:=CyclicClass(EE,QQ,dd);Append(~cls,cl);Append(~imageDegrees,dy);
   end for;
  else
   assert 3*dd le 7;
   K:=BaseRing(EE);LL<aa>:=ext<K|fac[1][1]>;
   EL:=EllipticCurve([LL!a:a in aInvariants(EE)]);PL:=EL![LL!PP[1],LL!PP[2]];
   TL:=EL![aa,0];QL:=8*PL+TL;assert 2*QL eq PL;
   cl,dy:=CyclicClass(EL,QL,3*dd);Append(~cls,cl);Append(~imageDegrees,dy);
  end if;
  mapCache[q]:=cls;
  Append(~localRecords,[*#localRecords+1,dd,split,imageDegrees,[Eltseq(a):a in cls]*]);
 end if;
 cls:=mapCache[q];
 if split then
  assert #cls eq 3;
  for m1 in [0..mul] do for m2 in [0..mul-m1] do
   m3:=mul-m1-m2;Include(~options,m1*cls[1]+m2*cls[2]+m3*cls[3]);
  end for;end for;
 else
  assert #cls eq 1;Include(~options,(mul div 3)*cls[1]);
 end if;
end procedure;
tested:=0;inertRejected:=0;classRejected:=0;survived:=0;classRecords:=[* *];survivors:=[* *];
for idx in [1..#targetdivisors] do
 TT:=targetdivisors[idx];ss,mm:=Support(TT);nc:=[i:i in [1..#ss]|ss[i] notin A3cs];
 if #nc eq 0 then continue;end if;
 tested+:=1;ee:=&+[mm[j]*Degree(ss[j]):j in nc];assert ee in [2..7];
 sums:={GY!0};parts:=[* *];
 for j in nc do
  opts:={GY!0};LocalOptions(ss[j],mm[j],~infoCache,~mapCache,~localRecords,~opts);
  Append(~parts,[*Degree(ss[j]),mm[j],#opts*]);
  sums:={a+b:a in sums,b in opts};
  if #sums eq 0 then break;end if;
 end for;
 if #sums eq 0 then
  inertRejected+:=1;
 else
  passes:=&or[qY(a) in complementsY[7-ee]:a in sums];
  Append(~classRecords,[*idx,ee,[Eltseq(a):a in sums],passes,parts*]);
  if passes then survived+:=1;Append(~survivors,[*idx,ee,parts*]);
  else classRejected+:=1;end if;
 end if;
 if tested mod 50 eq 0 then print "M230_D7_REVERSE_PROGRESS",tested,inertRejected,classRejected,survived,#mapCache,Cputime(t0);end if;
end for;
assert tested eq noncusp and tested eq inertRejected+classRejected+survived;
Write("data/M230_D7_reverse_arithmetic.txt",Sprint([*Invariants(GY),[Eltseq(r):r in Rows(RelationMatrix(GY))],[Eltseq(a):a in hgensY],dsY,[Eltseq(a):a in ccY],classRecords,localRecords,survivors*]) : Overwrite:=true);
printf "M230_D7_REVERSE_DONE target_classes=256 target_divisors=%o noncuspidal=%o inert_rejected=%o Yclass_rejected=%o remaining=%o complete=true time=%o\n",#targetdivisors,tested,inertRejected,classRejected,survived,Cputime(t0);
quit;
