/* Degree7 rational cuspidal-image filter on X1(2,24).
   All noncuspidal trace supports: single3..7,[3,3]+cusp1,[3,4].
   Exact integral-unit blocks lift Q; DEHMZProp4.16 identifies their
   subgroup with the full rational Jacobian after invariant check.
   Retain full cusp-complement classes, not merely their quotient cosets. */
SetColumns(0);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
Attach("fast_hecke.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("fast_diamond_hash_m2.m");
X:=MDX11(2,24,GF(5));C:=Curve(X);FF:=FunctionField(C);
b,c:=Explode([FF!f:f in X`_coordinates]);d:=b*(1-b-c);
E:=EllipticCurve([0,c,0,d,0]);Q:=E![b,b];fs:=[b,d,c^2-4*d];
for k in [1..12] do A:=k*Q;assert A[3] ne 0 and A[1] ne 0;Append(~fs,A[1]);end for;
cs:=Cusps(X);groups:=AssociativeArray();
for cusp in cs do
 sg:=[Valuation(f,cusp):f in fs];
 if IsDefined(groups,sg) then Append(~groups[sg],cusp);else groups[sg]:=[cusp];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));
blocks:=[&+[Divisor(cusp):cusp in groups[sg]]:sg in sigs];
degrees:=[Degree(o):o in blocks];
assert Sort(degrees) eq [1:i in [1..8]] cat [2:i in [1..6]] cat [4:i in [1..5]];
baseBlock:=Index(degrees,1);c0:=Support(blocks[baseBlock])[1];
print "N24_D7_SIGNATURES",#blocks,Sort(degrees),Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0);
A:=JH(48,[25]);assert Dimension(A) eq Genus(X);
assert &*[a:a in Invariants(G)|a ne 0] eq Integers()!Evaluate(FrobeniusPolynomial(A,5),1);
gens:=[m2(o-Degree(o)*Divisor(c0)):o in blocks];
Free:=FreeAbelianGroup(#blocks);h:=hom<Free->G|gens>;H:=Image(h);
assert Invariants(H) eq [2,4,4,120,240];QQ,qmap:=quo<G|H>;
print "N24_D7_GROUP",Invariants(G),Invariants(H),Invariants(QQ),Cputime(t0);
// These exact characteristic-zero witnesses are independently checked in
// verify_n24_reducible_witnesses.m. No degree6 exclusion is assumed.
R<V>:=PolynomialRing(Rationals());
fsQ:=[V^6-3*V^5+V^4/2+4*V^3-V^2-3*V/2+1/2,
 V^6+4*V^5+41*V^4/6+31*V^3/6+4*V^2/3-V/3-1/6];
hsQ:=[3-3*V-13*V^2+2*V^3+10*V^4-4*V^5,
 (-92-15*V+542*V^2+1074*V^3+744*V^4+216*V^5)/17];
Fp:=GF(5);S<T>:=PolynomialRing(Fp);known:={G!0};Exclude(~known,G!0);
for j in [1,2] do
 f:=S!fsQ[j];hh:=S!hsQ[j];assert IsIrreducible(f);
 K<aa>:=ext<Fp|f>;pt:=C(K)![Evaluate(hh,aa),aa,K!1];
 ps:=Places(pt);assert #ps eq 1 and Degree(ps[1]) eq 6;
 for a in [1,5,7,11] do
  xx:=DiamondOperator(X,a,ps[1]);cl:=m2(Divisor(xx)-6*Divisor(c0));assert cl in H;
  for bi in [1..#blocks] do
   if degrees[bi] eq 1 then Include(~known,cl+gens[bi]);end if;
  end for;
 end for;
end for;
assert #known eq 64;
print "N24_D7_KNOWN_REDUCIBLE",#known,Cputime(t0);
cuspclasses:=[m2(Divisor(cusp)-Degree(cusp)*Divisor(c0)):cusp in cs];
ds:=[Degree(cusp):cusp in cs];
procedure AddCuspClasses(remaining,start,current,~S)
 if remaining eq 0 then Include(~S,current);return;end if;
 for i in [start..#ds] do
  if ds[i] le remaining then AddCuspClasses(remaining-ds[i],i,current+cuspclasses[i],~S);end if;
 end for;
end procedure;
complements:=[* *];
for k in [0..4] do
 S:={G!0};Exclude(~S,G!0);AddCuspClasses(k,1,G!0,~S);aa:=AssociativeArray();
 for cl in S do
  ky:=-qmap(cl);
  if not IsDefined(aa,ky) then aa[ky]:={G!0};Exclude(~aa[ky],G!0);end if;
  Include(~aa[ky],cl);
 end for;
 Append(~complements,aa);
 print "N24_D7_COMPLEMENT",k,#Keys(aa),#S,Cputime(t0);
end for;
places:=[NoncuspidalPlaces(X,k):k in [1..7]];
assert [#xs:xs in places] eq [0,0,16,96,672,2712,11184];
liftclasses:={G!0};Exclude(~liftclasses,G!0);records:=[* *];tested:=0;passed:=0;
procedure TestClass(gc,lambda,indices,~liftclasses,~records,~tested,~passed)
 e:=&+lambda;k:=7-e;assert k in [0..4];tested+:=1;
 aa:=complements[k+1];ky:=qmap(gc);
 if IsDefined(aa,ky) then
  passed+:=1;
  for cc in aa[ky] do
   total:=gc+cc;assert total in H;Include(~liftclasses,total);
   Append(~records,[*lambda,indices,Eltseq(cc),Eltseq(total)*]);
  end for;
  print "N24_D7_PASS",lambda,indices,#aa[ky],passed,#liftclasses,Cputime(t0);
 end if;
end procedure;
repcounts:=[];
for dd in [3..7] do
 rs:=PlacesUpToDiamondM2HashDirect(X,places[dd]);passed0:=passed;tested0:=tested;
 assert #rs eq [4,32,168,684,2796][dd-2];Append(~repcounts,#rs);
 print "N24_D7_REPS",dd,#rs,Cputime(t0);
 for i in [1..#rs] do
  gc:=m2(Divisor(rs[i])-dd*Divisor(c0));
  TestClass(gc,[dd],[i],~liftclasses,~records,~tested,~passed);
  if i mod 100 eq 0 then print "N24_D7_PROGRESS",dd,i,#rs,passed-passed0,Cputime(t0);end if;
 end for;
 print "N24_D7_TYPE",[dd],tested-tested0,passed-passed0,Cputime(t0);
end for;
rs3:=PlacesUpToDiamondM2HashDirect(X,places[3]);
classes3:=[m2(Divisor(x)-3*Divisor(c0)):x in places[3]];
classes4:=[m2(Divisor(x)-4*Divisor(c0)):x in places[4]];
for ee in [3,4] do
 passed0:=passed;tested0:=tested;cls:=ee eq 3 select classes3 else classes4;
 for i in [1..#rs3] do
  gc:=m2(Divisor(rs3[i])-3*Divisor(c0));
  for j in [1..#cls] do
   TestClass(gc+cls[j],[3,ee],[i,j],~liftclasses,~records,~tested,~passed);
  end for;
 end for;
 assert tested-tested0 eq 4*[16,96][ee-2];
 print "N24_D7_TYPE",[3,ee],tested-tested0,passed-passed0,Cputime(t0);
end for;
assert tested eq 4+32+168+684+2796+64+384;
Kern:=Kernel(h);KB:=Matrix(Integers(),[Eltseq(Free!k):k in Generators(Kern)]);
L:=LatticeWithBasis(LLL(KB));vectors:=[];unknown:={gc:gc in liftclasses|gc notin known};
print "N24_D7_CLASS_STATUS",#liftclasses,#(liftclasses meet known),#unknown;
for gc in unknown do
 v:=Eltseq(gc@@h);w:=ClosestVector(L,Vector(Rationals(),v));
 short:=[v[i]-Integers()!w[i]:i in [1..#v]];
 assert h(Free!short) eq gc;Append(~vectors,short);
end for;
Write("data/N24_D7_cusp_lifts.txt",Sprint([*24,7,sigs,degrees,baseBlock,vectors*]) : Overwrite:=true);
Write("data/N24_D7_survivor_classes.txt",Sprint(records) : Overwrite:=true);
printf "N24_D7_DONE tested=%o retained=%o rational_classes=%o unresolved_classes=%o complete=true time=%o\n",tested,passed,#liftclasses,#unknown,Cputime(t0);
quit;
