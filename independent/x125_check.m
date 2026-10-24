// Independent check of the classification of septic points on X_1(25) (degree-7 draft, Prop. all-septimic-25),
// using the Sutherland--van Hoeij model (not the draft's code).
SetColumns(0); t0:=Cputime();
F2:=GF(2);
// ---- model over F_2 -------------------------------------------------------------
Fu<u>:=RationalFunctionField(F2); P<V>:=PolynomialRing(Fu);
Xpol:=u*V^5 + (u^4 - 2*u^3 - u^2 + 2*u + 1)*V^4 - (2*u^6 - 2*u^4 + 4*u^3 + 2*u^2 - 2)*V^3
     + (u^8 + u^7 - 2*u^6 + u^5 - u^4 - u^3 - 2*u^2 - u + 1)*V^2 + (u^8 + u^7 + 2*u^6 + u^5 - 2*u^4 + u^3 - u^2)*V + u^6;
assert IsIrreducible(Xpol);
F<v>:=FunctionField(Xpol);
y:=(u-v-1)/u;
x:=-(v+1)*(u^4 + u^3*v + u^3 - u*v^2 - v)/(u^2*(u^3 + u^2*v + u^2 - u*v - u - v^2 - v));
r:=(x^2*y - x*y + y - 1)/(x^2*y - x); s:=(x*y - y + 1)/(x*y);
b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); P0:=E![0,0];
g:=Genus(F); j:=jInvariant(E);
printf "genus %o, deg j %o\n", g, Degree(j);
assert g eq 12 and Degree(j) eq 300;
Q5:=5*P0; assert Q5 ne E!0; assert 5*Q5 eq E!0;
print "base point has exact order 25 over the function field";
// ---- cusps, places, class group ----------------------------------------------------
cusps:=[Q: Q in Support(Divisor(j)) | Valuation(j,Q) lt 0];
printf "cusp degrees %o\n", Sort([Degree(Q): Q in cusps]);
nonc:=[[Q: Q in Places(F,d) | Valuation(j,Q) ge 0]: d in [1..7]];
printf "noncuspidal place counts %o (%o s)\n", [#L: L in nonc], Cputime(t0);
Cl,m:=ClassGroup(F); printf "class group %o (%o s)\n", Invariants(Cl), Cputime(t0);
c1:=[Q: Q in cusps | Degree(Q) eq 1]; c0:=c1[1];
cl:=func<D | (D - Degree(D)*c0) @@ m>;
// all effective cuspidal divisors of degree 7 and their classes
degs:=[Degree(Q): Q in cusps];
procedure rec(i, rem, cur, ~L) if rem eq 0 then Append(~L,cur); return; end if; if i gt #cusps then return; end if;
  for k in [0..(rem div degs[i])] do rec(i+1, rem-k*degs[i], cur+k*cusps[i], ~L); end for; end procedure;
L7:=[]; rec(1,7,DivisorGroup(F)!0,~L7);
cusp7:={cl(D): D in L7}; printf "effective cuspidal degree-7 divisors %o, classes %o\n", #L7, #cusp7;
L2:=[]; rec(1,2,DivisorGroup(F)!0,~L2);
bad5:=[<Q,C>: Q in nonc[5], C in L2 | cl(1*Q+C) notin cusp7];
printf "degree-5 place + cuspidal degree-2: %o divisors, %o with class outside the cuspidal classes\n", #nonc[5]*#L2, #bad5;
// ---- reductions of the diamond translates of the known sextic and septic points ---------
function tate(E,Q) // Tate normal form parameters (b,c) of (E,Q), Q of order >= 4
  a1,a2,a3,a4,a6:=Explode(aInvariants(E)); x0:=Q[1]; y0:=Q[2];
  // translate Q to (0,0)
  E1:=Transformation(E,[x0,0,y0,1]); a1,a2,a3,a4,a6:=Explode(aInvariants(E1));
  assert a6 eq 0; E2:=Transformation(E1,[0,a4/a3,0,1]); a1,a2,a3,a4,a6:=Explode(aInvariants(E2));
  assert a4 eq 0 and a6 eq 0; bb:=-a2^3/a3^2; cc:=1-a1*a2/a3; return bb,cc;
end function;
function key(bb,cc) // Frobenius orbit of (b,c), as a set of pairs in a fixed field
  G:=GF(2,LCM(Degree(MinimalPolynomial(bb,F2)),Degree(MinimalPolynomial(cc,F2))));
  return {<MinimalPolynomial(bb,F2),MinimalPolynomial(cc,F2)>}, G; end function;
Qx<X>:=PolynomialRing(Rationals());
function reduce_point(f, ypol_or_func, deg)
  return 0; end function;
// sextic point: x^3-x^2+1=0, y^2+(x^2-2x-1)y+x=0 ; septic point: a^7+a^6-a^5-a^4+a^2+a-1=0, y=a^2+a
K3<w>:=NumberField(X^3-X^2+1); K3y<Y>:=PolynomialRing(K3); K6<z>:=NumberField(Y^2+(w^2-2*w-1)*Y+w);
K6a:=AbsoluteField(K6); xs:=K6a!w; ys:=K6a!z;
K7<a>:=NumberField(X^7+X^6-X^5-X^4+X^2+X-1); x7:=a; y7:=a^2+a;
function points_over_F2(K, xx, yy)
  rr:=(xx^2*yy - xx*yy + yy - 1)/(xx^2*yy - xx); ss:=(xx*yy - yy + 1)/(xx*yy);
  bb:=rr*ss*(rr-1); cc:=ss*(rr-1); EK:=EllipticCurve([1-cc,-bb,-bb,0,0]); PK:=EK![0,0];
  assert Order(PK) eq 25;
  O:=MaximalOrder(K); pr:=Factorization(2*O);
  out:=[];
  for pp in pr do kk,red:=ResidueClassField(pp[1]);
    Ek:=EllipticCurve([red(1-cc),red(-bb),red(-bb),0,0]); Pk:=Ek![0,0]; assert Order(Pk) eq 25;
    for aa in [1,2,3,4,6,7,8,9,11,12] do b1,c1:=tate(Ek,aa*Pk); Append(~out,<Degree(kk),b1,c1>); end for;
  end for; return out;
end function;
function place_key(Q) k,mk:=ResidueClassField(Q); return <Degree(Q), MinimalPolynomial(mk(b),F2), MinimalPolynomial(mk(c),F2), mk(b), mk(c)>; end function;
function samepoint(k, bq, cq, kp, b1, c1) // is there a field iso k -> kp sending (bq,cq) to (b1,c1)?
  if #k ne #kp then return false; end if;
  f:=DefiningPolynomial(k); PP:=PolynomialRing(kp); 
  for rt in Roots(PP!Eltseq(f)) do beta:=rt[1];
    phi:=func<z | &+[ kp!(Eltseq(z)[i]) * beta^(i-1) : i in [1..#Eltseq(z)]]>;
    if phi(bq) eq b1 and phi(cq) eq c1 then return true; end if; end for;
  return false; end function;
function matchplace(L, tri) d,b1,c1:=Explode(tri); hits:=[];
  for Q in L do k,mk:=ResidueClassField(Q); if Degree(Q) ne d then continue; end if;
    if samepoint(k, mk(b), mk(c), Parent(b1), b1, c1) then Append(~hits,Q); end if; end for;
  return hits; end function;
sext:=points_over_F2(K6a, xs, ys); sep:=points_over_F2(K7, x7, y7);
printf "residue degrees at 2: sextic %o, septic %o\n", {t[1]: t in sext}, {t[1]: t in sep};
cov6:={}; for t in sext do if t[1] eq 6 then h:=matchplace(nonc[6],t); assert #h eq 1; Include(~cov6,h[1]); end if; end for;
cov7:={}; for t in sep do if t[1] eq 7 then h:=matchplace(nonc[7],t); assert #h eq 1; Include(~cov7,h[1]); end if; end for;
printf "sextic places covered by reductions of translates: %o of %o\n", #cov6, #nonc[6];
printf "septic places covered by reductions of translates: %o of %o\n", #cov7, #nonc[7];
printf "X125_CHECK_DONE (%o s)\n", Cputime(t0);
quit;
