/* Exact finite-class necessary conditions for closed septimic X1(25)
   points. No degree-seven function exists (DvH Theorem3). Atp2 all
   closed cusps liftQ, since2 generates (Z/25)^*. J25(Q) is odd.
   Therefore a noncuspidal degree7 divisor cannot share its Q class
   with an effective cuspidal divisor. */
SetColumns(0);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
X:=MDX1(25,GF(2));C:=Curve(X);cs:=Cusps(X);ds:=[Degree(c):c in cs];
assert Genus(X) eq 12 and Sort(ds) eq [1:i in [1..10]] cat [4,4,10];
c0:=[c:c in cs|Degree(c) eq 1][1];
print "C25_P2_SETUP",Genus(X),ds,Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
f:=FrobeniusPolynomial(JOne(25),2);assert Degree(f) eq 24;
assert &*[a:a in Invariants(G)|a ne 0] eq Integers()!Evaluate(f,1);
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
HQ:=sub<G|cc>;assert Invariants(HQ) eq [227555];
print "C25_P2_CLASSGROUP",Invariants(G),Invariants(HQ),Cputime(t0);
procedure AddEffective(remaining,start,current,~out,~number)
 if remaining eq 0 then Include(~out,current);number+:=1;return;end if;
 for i in [start..#cs] do
  if ds[i] le remaining then AddEffective(remaining-ds[i],i,current+cc[i],~out,~number);end if;
 end for;
end procedure;
complements:=AssociativeArray();
for k in [0..2] do
 out:={G|};number:=0;AddEffective(k,1,G!0,~out,~number);complements[k]:=out;
 print "C25_P2_COMPLEMENTS",k,number,#out,Cputime(t0);
end for;
qeffective:={G|};number:=0;AddEffective(7,1,G!0,~qeffective,~number);
print "C25_P2_EFFECTIVE_CUSPIDAL",number,#qeffective,Cputime(t0);
places:=[NoncuspidalPlaces(X,d):d in [1..7]];
assert [#places[d]:d in [1..4]] eq [0,0,0,0];
print "C25_P2_NONCUSP",[#ss:ss in places],Cputime(t0);
traces:=[];closed:=[];
for n in [1..7] do
 sn:=-n*Coefficient(f,24-n);
 for j in [1..n-1] do sn-:=Coefficient(f,24-j)*traces[n-j];end for;
 Append(~traces,Integers()!sn);rational:=2^n+1-Integers()!sn;
 proper:=[d:d in Divisors(n)|d lt n];
 count:=rational;for d in proper do count-:=d*closed[d];end for;
 assert count mod n eq 0;Append(~closed,count div n);
 assert closed[n] eq #places[n]+#[c:c in cs|Degree(c) eq n];
end for;
print "C25_P2_COUNTS_FROBENIUS_VERIFIED",closed;
candidates:={G|};records:=[* *];
for d in [5..7] do
 passes:=0;effective:=0;
 for idx in [1..#places[d]] do
  x:=places[d][idx];a:=m2(Divisor(x)-d*Divisor(c0));
  for b in complements[7-d] do
   gc:=a+b;if gc notin HQ then continue;end if;passes+:=1;
   if gc in qeffective then effective+:=1;continue;end if;
   Include(~candidates,gc);Append(~records,[*d,idx,Eltseq(a),Eltseq(b),Eltseq(gc)*]);
  end for;
 end for;
 print "C25_P2_TYPE",d,#places[d],passes,effective,#candidates,Cputime(t0);
end for;
print "C25_P2_CANDIDATE_CLASSES",[Eltseq(a):a in candidates];
print "C25_P2_CANDIDATE_RECORDS",records;
printf "C25_P2_CLASS_FILTER_DONE candidates=%o records=%o complete=true time=%o\n",#candidates,#records,Cputime(t0);
/* Retain the exact group/maps for a continuation comparing known points. */
