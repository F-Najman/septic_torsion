/* Complete first-support diamond-orbit reduction of the degree7 class test.
   Load after probe_36_class5_v2.m. HQ, Q-effective cusp classes, and
   the known global orbit union are diamond-stable. Every finite
   divisor has a diamond translate considered below. */
assert Invariants(HQ) eq [12,252,7812];
assert [#ss:ss in places][1..6] eq [0,0,48,144,660,2856];
trace36:=[];closed36:=[];
for n36 in [1..7] do
 sn36:=-n36*Coefficient(fp,34-n36);
 for j36 in [1..n36-1] do sn36-:=Coefficient(fp,34-j36)*trace36[n36-j36];end for;
 Append(~trace36,Integers()!sn36);count36:=5^n36+1-Integers()!sn36;
 for d36 in Divisors(n36) do
  if d36 lt n36 then count36-:=d36*closed36[d36];end if;
 end for;
 assert count36 mod n36 eq 0;Append(~closed36,count36 div n36);
 assert closed36[n36] eq #places[n36]+#[c:c in cs|Degree(c) eq n36];
end for;
print "C36_FROBENIUS_COUNTS_VERIFIED",[#ss:ss in places],closed36,Cputime(t0);
procedure C36OrbitEffective(weights,classes,remaining,start,current,~out,~number)
 if remaining eq 0 then Include(~out,current);number+:=1;return;end if;
 for i in [start..#weights] do
  if weights[i] le remaining then
   C36OrbitEffective(weights,classes,remaining-weights[i],i,current+classes[i],~out,~number);
  end if;
 end for;
end procedure;
qeffective:={G|};number:=0;
C36OrbitEffective([Degree(o):o in os],qc,7,1,G!0,~qeffective,~number);
assert number eq 19159 and #qeffective eq number;
print "C36_DEGREE7_Q_EFFECTIVE",number,#qeffective,Cputime(t0);
allcusps:={G|};number:=0;
C36OrbitEffective([Degree(c):c in cs],cc,7,1,G!0,~allcusps,~number);
assert number eq 188865;
candidates:={gc:gc in allcusps|gc in HQ and gc notin qeffective};
print "C36_DEGREE7_ALL_CUSPIDAL",number,#allcusps,"new_candidates",#candidates,Cputime(t0);
Qbar,qbar:=quo<G|HQ>;complements:=AssociativeArray();
for k in [0..4] do
 out:={G|};number:=0;
 C36OrbitEffective([Degree(c):c in cs],cc,k,1,G!0,~out,~number);
 bycoset:=AssociativeArray();
 for a in out do
  key:=qbar(a);
  if IsDefined(bycoset,key) then Include(~bycoset[key],a);
  else bycoset[key]:={a};end if;
 end for;
 complements[k]:=bycoset;
 print "C36_DEGREE7_COMPLEMENTS",k,number,#out,#Keys(bycoset),Cputime(t0);
end for;
Attach("fast_diamond_hash_m2.m");
Attach("fast_diamond_hash_cyclic_bc.m");
reps36:=[];
for degree36 in [1..7] do
 if IsEmpty(places[degree36]) then rr36:=[];
 else rr36:=PlacesUpToDiamondCyclicHashBC(X,places[degree36]);end if;
 Append(~reps36,rr36);
 print "C36_ORBIT_REPRESENTATIVES",degree36,#places[degree36],#rr36,Cputime(t0);
end for;
assert #reps36[3] eq 8 and #reps36[7] eq 1788;
repclasses36:=[];
for low36 in [1..5] do
 Append(~repclasses36,[m2(Divisor(x)-low36*Divisor(c0)):x in reps36[low36]]);
end for;
print "C36_ORBIT_LOW_CLASSES_READY",[#s:s in repclasses36],Cputime(t0);
