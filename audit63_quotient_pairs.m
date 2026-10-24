SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
t0:=Cputime();
function E9at(t)
 b:=t^2*(t-1)*(t^2-t+1);c:=t^2*(t-1);K:=Parent(t);
 return EllipticCurve([1-c,-b,-b,K!0,K!0]);
end function;
function NormalizePair(E,R,S)
 assert R[3] eq 1 and S[3] eq 1;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x0:=R[1];y0:=R[2];
 aa3:=2*y0+a3+a1*x0;aa2:=3*x0+a2;aa4:=3*x0^2+2*x0*a2+a4-a1*y0;
 slope:=aa4/aa3;atwo:=aa2-slope^2-a1*slope;scale:=aa3/atwo;
 b:=-atwo^3/aa3^2;c:=-((2*slope+a1)*atwo-aa3)/aa3;
 t:=c^2/(b-c);x:=(S[1]-x0)/scale^2;y:=(S[2]-y0-slope*(S[1]-x0))/scale^3;
 En:=E9at(t);assert aInvariants(En) eq [1-c,-b,-b,Parent(t)!0,Parent(t)!0];
 Q:=En![x,y];assert 7*Q eq En!0 and Q ne En!0;
 return t,x;
end function;
F<z>:=GF(128);pairs:=[];
for tt in F do
 if tt in [F!0,F!1] then continue;end if;
 try E:=E9at(tt);catch e;continue;end try;
 
 for rr in Roots(DivisionPolynomial(E,7)) do
  assert #Points(E,rr[1]) eq 2 and #E eq 126;
  Append(~pairs,<tt,rr[1]>);
 end for;
end for;
assert #pairs eq 189;
function PairKey(t,x)
 return Minimum([Eltseq(t^(2^j),GF(2)) cat Eltseq(x^(2^j),GF(2)):j in [0..6]]);
end function;
lookup:=AssociativeArray();for r in pairs do lookup[PairKey(r[1],r[2])]:=r;end for;
assert #Keys(lookup) eq 27;
left:=Keys(lookup);reps:=[];
while #left gt 0 do
 key:=Minimum(Setseq(left));Append(~reps,lookup[key]);todo:=[key];orb:={key};
 while #todo gt 0 do
  k:=todo[1];Remove(~todo,1);r:=lookup[k];E:=E9at(r[1]);P:=E![0,0];Q:=Points(E,r[2])[1];
  for a in [2,5] do
   ta,xa:=NormalizePair(E,a*P,a*Q);ka:=PairKey(ta,xa);assert IsDefined(lookup,ka);
   if ka notin orb then Include(~orb,ka);Append(~todo,ka);end if;
  end for;
 end while;
 assert #orb eq 9;left diff:=orb;
end while;
assert #reps eq 3;
printf "Y63_NC7_EXHAUSTIVE geometric_pairs=%o closed_places=%o diamond_orbits=%o time=%o\n",#pairs,#Keys(lookup),#reps,Cputime(t0);
for i in [1..#reps] do printf "Y63_NC7_REP index=%o fieldpoly=%o t=%o x7=%o\n",i,Eltseq(DefiningPolynomial(F,GF(2))),Eltseq(reps[i][1],GF(2)),Eltseq(reps[i][2],GF(2));end for;
print "Y63_ALL_QUOTIENT_NC7_PAIRS_VERIFIED";quit;
