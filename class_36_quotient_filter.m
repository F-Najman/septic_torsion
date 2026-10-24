/* Complete target classification after a necessary genus7 class test.
   Inputs: probe36, prepare36, quotient36. Every first support is
   reduced by diamonds; all second supports and cusps remain. */
assert [#s:s in reps36] eq [0,0,8,24,110,480,1788];
assert [#s:s in repclasses36] eq [0,0,8,24,110];
assert Invariants(HQ) eq [12,252,7812] and Invariants(HY36) eq [2,42,126];
xcache36:=AssociativeArray();
if assigned precomputed6classes36 then
 assert #precomputed6classes36 le #reps36[6];
 for ii36 in [1..#precomputed6classes36] do
  xcache36[reps36[6][ii36]]:=precomputed6classes36[ii36];
 end for;
end if;
smallclasses36:=[[m2(Divisor(x)-low36*Divisor(c0)):x in places[low36]]:low36 in [3,4]];
candidates:={gc:gc in allcusps|gc in HQ and gc notin qeffective};assert IsEmpty(candidates);
procedure C36QuotientFineTest(a,degree,HQ,qbar,complements,qeffective,~candidates,~passes,~effective)
 key:=-qbar(a);tbl:=complements[7-degree];
 if not IsDefined(tbl,key) then return;end if;
 for b in tbl[key] do
  gc:=a+b;assert gc in HQ;passes+:=1;
  if gc in qeffective then effective+:=1;else Include(~candidates,gc);end if;
 end for;
end procedure;
for low36 in [3..5] do
 passes36:=0;effective36:=0;
 for a36 in repclasses36[low36] do
  C36QuotientFineTest(a36,low36,HQ,qbar,complements,qeffective,~candidates,~passes36,~effective36);
 end for;
 print "C36_QFILTER_LOW_TYPE",[low36],#repclasses36[low36],passes36,effective36,#candidates,Cputime(t0);
end for;
for second36 in [3,4] do
 passes36:=0;effective36:=0;total36:=0;
 for a36 in repclasses36[3] do
  for b36 in smallclasses36[second36-2] do
   total36+:=1;
   C36QuotientFineTest(a36+b36,3+second36,HQ,qbar,complements,qeffective,~candidates,~passes36,~effective36);
  end for;
 end for;
 assert total36 eq #reps36[3]*#places[second36];
 print "C36_QFILTER_LOW_TYPE",[3,second36],total36,passes36,effective36,#candidates,Cputime(t0);
end for;
highrecords36:=[* *];
journal36:=Open("logs/x1_36/quotient_filter.log","w");
fprintf journal36,"C36_QFILTER_INPUT_COUNTS %o representatives %o\n",[#s:s in places],[#s:s in reps36];
Flush(journal36);
for high36 in [6,7] do
 ypasses36:=0;finepasses36:=0;effective36:=0;
 for ihigh36 in [1..#reps36[high36]] do
  source36:=reps36[high36][ihigh36];
  yyclass36:=C36YClass(X,Y36,source36,my36,c0y36);
  oky36:=-qy36(yyclass36) in yc36[7-high36];
  xxrecord36:=[];
  if oky36 then
   ypasses36+:=1;
   if IsDefined(xcache36,source36) then xxclass36:=xcache36[source36];
   else xxclass36:=m2(Divisor(source36)-high36*Divisor(c0));xcache36[source36]:=xxclass36;end if;
   xxrecord36:=Eltseq(xxclass36);
   C36QuotientFineTest(xxclass36,high36,HQ,qbar,complements,qeffective,~candidates,~finepasses36,~effective36);
  end if;
  Append(~highrecords36,[*high36,ihigh36,Eltseq(yyclass36),oky36,xxrecord36*]);
  fprintf journal36,"C36_QFILTER_CASE degree=%o index=%o ypass=%o yc=%o xc=%o candidates=%o time=%o\n",
   high36,ihigh36,oky36,Eltseq(yyclass36),xxrecord36,#candidates,Cputime(t0);
  Flush(journal36);
  if ihigh36 mod 20 eq 0 or ihigh36 eq #reps36[high36] then
   print "C36_QFILTER_PROGRESS",high36,ihigh36,#reps36[high36],ypasses36,#candidates,Cputime(t0);
  end if;
 end for;
 print "C36_QFILTER_HIGH_TYPE",[high36],#reps36[high36],ypasses36,finepasses36,effective36,#candidates,Cputime(t0);
end for;
assert #highrecords36 eq 480+1788;
fprintf journal36,"C36_QFILTER_COMPLETE records=%o candidates=%o complete=true\n",#highrecords36,#candidates;
Flush(journal36);delete journal36;
print "C36_QUOTIENT_ORBIT_FILTER_DONE records",#highrecords36,"candidates",#candidates,"complete=true",Cputime(t0);
