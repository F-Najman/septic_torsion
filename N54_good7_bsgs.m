/* Complete X1(54) pure-degree7 image orbits on X1(27)/F5.
   Exact Hasse BSGS counts and explicit prime-degree coordinate ideals.
   No ordinary #E or Order call occurs in the production loop. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("prime_degree_place.m");
Attach("hasse_bsgs_count.m");
X:=MDX1(27,GF(5));FF:=FunctionField(Curve(X));K:=GF(5^7);
assert (#K-1) mod 3 ne 0;
// Exercise the swapped coordinate branch without creating a subfield.
f,g:=MDTwoGeneratorsPrimeDegree(K!2,K.1);
assert Evaluate(f,[K!2,K.1]) eq 0 and Evaluate(g,[K!2,K.1]) eq 0;
seen:={K|};out:=[];records:=[* *];weighted:=0;tot:=0;jcount:=0;curves:=0;
bsgsTime:=0;mapTime:=0;maxTrials:=0;rootHistogram:=AssociativeArray();
rootHistogram[1]:=0;rootHistogram[3]:=0;
for j in K do
 if j in seen then continue;end if;
 for i in [0..6] do Include(~seen,j^(5^i));end for;
 if Degree(MinimalPolynomial(j,GF(5))) ne 7 then continue;end if;
 jcount+:=1;es:=Twists(EllipticCurveWithjInvariant(j));assert #es eq 2;
 for ie in [1..2] do
  E:=es[ie];assert BaseRing(E) cmpeq K;
  t:=Cputime();n,trials:=EllipticCountHasseBSGS(E);bsgsTime+:=Cputime(t);
  maxTrials:=Max(maxTrials,trials);curves+:=1;
  if n mod 54 ne 0 then continue;end if;
  a:=Valuation(n,3);assert a ge 3;found:=false;
  for trial in [1..100] do
   Q:=(n div 3^a)*Random(E);assert (3^a)*Q eq E!0;
   if (3^(a-1))*Q ne E!0 then found:=true;break;end if;
  end for;
  assert found;P:=3^(a-3)*Q;assert 27*P eq E!0 and 9*P ne E!0;
  xy:=TateNormalForm_xy(E,P);
  t:=Cputime();p:=MDPlaceCommonZerosPrimeDegree([FF.1,FF.2],xy);mapTime+:=Cputime(t);
  assert Degree(p) eq 7;
  if tot eq 0 then
   ref:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=P>);
   assert p eq ref;print "N54_BSGS_FIRST_INVERSE_REFERENCE_PASSED",Cputime(t0);
  end if;
  b2,b4,b6,b8:=Explode(bInvariants(E));R<x>:=PolynomialRing(K);
  roots:=Roots(4*x^3+b2*x^2+2*b4*x+b6);r:=#roots;assert r in [1,3];
  Append(~out,p);weighted+:=9*r;tot+:=1;rootHistogram[r]+:=1;
  Append(~records,[*[[Integers()!c:c in Eltseq(z,GF(5))]:z in xy],n,r,
                    [Integers()!c:c in Eltseq(j,GF(5))]*]);
  if tot mod 50 eq 0 then print "N54_BSGS_GOOD7_PROGRESS",tot,weighted,jcount,bsgsTime,mapTime,Cputime(t0);end if;
 end for;
 if jcount mod 500 eq 0 then print "N54_BSGS_J_PROGRESS",jcount,curves,tot,weighted,Cputime(t0);end if;
end for;
assert #seen eq #K and jcount eq (5^7-5) div 7 and curves eq 2*jcount;
assert #Seqset(out) eq #out and #out eq #records;
assert weighted eq 11412 and weighted eq 9*(rootHistogram[1]+3*rootHistogram[3]);
fieldpolynomial:=[Integers()!c:c in Coefficients(DefiningPolynomial(K,GF(5)))];
Write("data/N54_D7_good7_pointdata.txt",Sprint([*fieldpolynomial,records*]) : Overwrite:=true);
printf "N54_BSGS_GOOD7_DONE representatives=%o weighted_source_places=%o j_orbits=%o curves=%o roots1=%o roots3=%o max_count_trials=%o complete=true time=%o\n",#out,weighted,jcount,curves,rootHistogram[1],rootHistogram[3],maxTrials,Cputime(t0);
print "N54_BSGS_TIMING count",bsgsTime,"inverse",mapTime;
quit;
