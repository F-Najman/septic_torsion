// Shared helpers for the independent degree-7 checks.
function tate(E,Q) // Tate normal form parameters (b,c) of (E,Q), Q of order >= 4, by direct substitution
  k:=BaseRing(E); a1,a2,a3,a4,a6:=Explode(aInvariants(E)); x0:=Q[1]; y0:=Q[2];
  R<X,Y>:=PolynomialRing(k,2); f:=Y^2+a1*X*Y+a3*Y-(X^3+a2*X^2+a4*X+a6);
  g:=Evaluate(f,[X+x0,Y+y0]); assert MonomialCoefficient(g,R!1) eq 0;
  A3:=MonomialCoefficient(g,Y); A4:=-MonomialCoefficient(g,X); assert A3 ne 0;
  g:=Evaluate(g,[X,Y+(A4/A3)*X]); assert MonomialCoefficient(g,X) eq 0;
  A1:=MonomialCoefficient(g,X*Y); A2:=-MonomialCoefficient(g,X^2); A3:=MonomialCoefficient(g,Y);
  assert MonomialCoefficient(g,Y^2) eq 1 and MonomialCoefficient(g,X^3) eq -1;
  return -A2^3/A3^2, 1-A1*A2/A3;
end function;
function ptkey(bb, cc) // canonical identity of the closed point with coordinates (bb,cc), independent of the field presentation
  k:=Parent(bb); F0:=PrimeField(k); n:=Degree(k,F0); g:=bb; h:=cc; flag:=true;
  if Degree(MinimalPolynomial(g,F0)) lt n then g:=cc; h:=bb; flag:=false; end if;
  assert Degree(MinimalPolynomial(g,F0)) eq n;
  M:=Matrix(F0,[Eltseq(g^i,F0): i in [0..n-1]]); sol:=Solution(M,Vector(F0,Eltseq(h,F0)));
  return <Eltseq(MinimalPolynomial(g,F0)), Eltseq(sol), flag>;
end function;
function samepoint(k, bq, cq, kp, b1, c1) // do (bq,cq) and (b1,c1) define the same closed point?
  if #k ne #kp then return false; end if; return ptkey(bq,cq) eq ptkey(b1,c1); end function;
function matchplace(L, tri, b, c) d,b1,c1:=Explode(tri); hits:=[];
  for Q in L do if Degree(Q) ne d then continue; end if; k,mk:=ResidueClassField(Q);
    if samepoint(k, mk(b), mk(c), Parent(b1), b1, c1) then Append(~hits,Q); end if; end for;
  return hits; end function;
procedure effdivs(cusps, i, rem, cur, ~L) // all effective divisors of degree rem supported on cusps[i..]
  if rem eq 0 then Append(~L,cur); return; end if; if i gt #cusps then return; end if;
  for k in [0..(rem div Degree(cusps[i]))] do effdivs(cusps, i+1, rem-k*Degree(cusps[i]), cur+k*cusps[i], ~L); end for;
end procedure;
// reductions at primes above p of the diamond translates of a known point given in raw (x,y) coordinates
function known_reductions(K, xx, yy, N, p)
  rr:=(xx^2*yy - xx*yy + yy - 1)/(xx^2*yy - xx); ss:=(xx*yy - yy + 1)/(xx*yy);
  bb:=rr*ss*(rr-1); cc:=ss*(rr-1); EK:=EllipticCurve([1-cc,-bb,-bb,0,0]); assert Order(EK![0,0]) eq N;
  O:=MaximalOrder(K); out:=[* *];
  reps:=[a: a in [1..N div 2] | Gcd(a,N) eq 1];
  for pp in Factorization(p*O) do assert pp[2] eq 1; kk,red:=ResidueClassField(pp[1]);
    Ek:=EllipticCurve([red(1-cc),red(-bb),red(-bb),0,0]); Pk:=Ek![0,0]; assert Order(Pk) eq N;
    Append(~out, [ <Degree(kk), b1, c1> where b1,c1:=tate(Ek, a*Pk) : a in reps ]);
  end for; return out; // one list of translates per prime above p
end function;
