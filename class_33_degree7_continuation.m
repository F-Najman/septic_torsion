/* Load in the completed probe_33_class2.m session. All effective F2
   degree7 divisors are considered, including the all-cuspidal ones. */
assert Invariants(G) eq [15,42373650,0] and Invariants(HQ) eq [5,42373650];
assert [#ss:ss in places] eq [0,0,0,0,2,20,30];
assert 3*#HQ eq &*[n:n in Invariants(G)|n ne 0];
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
procedure C33Effective(weights,classes,remaining,start,current,~out,~number)
 if remaining eq 0 then Include(~out,current);number+:=1;return;end if;
 for i in [start..#weights] do
  if weights[i] le remaining then
   C33Effective(weights,classes,remaining-weights[i],i,current+classes[i],~out,~number);
  end if;
 end for;
end procedure;
qeffective:={G|};number:=0;
C33Effective([Degree(o):o in os],qc,7,1,G!0,~qeffective,~number);
assert number eq 25100 and #qeffective eq 25100;
print "C33_DEGREE7_Q_EFFECTIVE",number,#qeffective,Cputime(t0);
allcusps:={G|};number:=0;
C33Effective([Degree(c):c in cs],cc,7,1,G!0,~allcusps,~number);
assert number eq 25220 and #allcusps eq 25220;
print "C33_DEGREE7_ALL_CUSPIDAL",number,#allcusps,Cputime(t0);
candidates:={gc:gc in allcusps|gc in HQ and gc notin qeffective};
assert #candidates eq 0;
print "C33_DEGREE7_ALL_CUSPIDAL_CANDIDATES",#candidates,Cputime(t0);
complements:=AssociativeArray();
for k in [0..2] do
 out:={G|};number:=0;
 C33Effective([Degree(c):c in cs],cc,k,1,G!0,~out,~number);complements[k]:=out;
 print "C33_DEGREE7_COMPLEMENTS",k,number,#out,Cputime(t0);
end for;
for d in [5..7] do
 total:=0;passes:=0;effective:=0;
 for x in places[d] do
  a:=m2(Divisor(x)-d*Divisor(c0));
  for b in complements[7-d] do
   total+:=1;gc:=a+b;if gc notin HQ then continue;end if;passes+:=1;
   if gc in qeffective then effective+:=1;continue;end if;
   Include(~candidates,gc);
  end for;
 end for;
 assert total eq [120,200,30][d-4] and passes eq [0,0,10][d-4] and effective eq 0;
 print "C33_DEGREE7_NONCUSP_TYPE",d,#places[d],total,passes,effective,#candidates,Cputime(t0);
end for;
print "C33_DEGREE7_CANDIDATES",[Eltseq(gc):gc in candidates];
/* The archived septic point has irreducible reduction at2. */
L:=GF(2^7);RL<z>:=PolynomialRing(L);
aa:=Roots(z^7+z^5+z^3+z+1)[1][1];yy:=aa^6+aa^4+aa^3+aa^2+aa;
pp:=Places(C(L)![aa,yy]);assert #pp eq 1 and Degree(pp[1]) eq 7 and not IsCusp(X,pp[1]);
known:=Seqset(DiamondOrbit(X,pp[1]));
knownclasses:={m2(Divisor(x)-7*Divisor(c0)):x in known};
print "C33_DEGREE7_KNOWN_ORBIT",#known,#knownclasses,"remaining_unknown",#(candidates diff knownclasses),Cputime(t0);
assert #known eq 10 and #knownclasses eq 10 and knownclasses subset candidates;
unknown:=candidates diff knownclasses;
assert #candidates eq 10 and #unknown eq 0;
printf "C33_DEGREE7_CLASS_COMPARISON_DONE candidates=%o known_classes=%o unknown_classes=%o complete=true time=%o\n",#candidates,#knownclasses,#unknown,Cputime(t0);
/* Retain unknown classes if a cross-prime RR continuation is needed. */
