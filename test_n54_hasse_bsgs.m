SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
Attach("hasse_bsgs_count.m");
K:=GF(5^7);seen:={K|};jcount:=0;curveCount:=0;bsgsTime:=0;refTime:=0;trials:=[];
for j in K do
 if j in seen then continue;end if;
 for i in [0..6] do Include(~seen,j^(5^i));end for;
 if Degree(MinimalPolynomial(j,GF(5))) ne 7 then continue;end if;
 jcount+:=1;
 for E in Twists(EllipticCurveWithjInvariant(j)) do
  t:=Cputime();n,r:=EllipticCountHasseBSGS(E);bsgsTime+:=Cputime(t);
  Append(~trials,r);curveCount+:=1;
  if curveCount le 100 then
   t:=Cputime();ref:=#E;refTime+:=Cputime(t);assert n eq ref;
  end if;
  if curveCount mod 20 eq 0 then print "N54_BSGS_TEST_PROGRESS",curveCount,bsgsTime,refTime,Max(trials),Cputime(t0);end if;
 end for;
 if jcount ge 100 then break;end if;
end for;
assert curveCount eq 200 and #trials eq 200;
print "N54_BSGS_TEST_DONE curves=200 exact_reference_comparisons=100 complete=true",bsgsTime,refTime,Max(trials),Cputime(t0);quit;
