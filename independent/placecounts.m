// Closed noncuspidal places of degree d<=D on X_1(n) and X_1(2,n) over F_p, by enumerating
// all elliptic curves over F_{p^d} (j-invariants and all twists) and counting level structures
// modulo automorphisms (free for level >= 5).  Independent of any model of the modular curves.
SetColumns(0);
if Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if Type(D) eq MonStgElt then D:=StringToInteger(D); end if;
Ns:=[StringToInteger(s): s in Split(Ns,",")];
Ms:=assigned M2s select [StringToInteger(s): s in Split(M2s,",")] else [];
function exact(a,b,n) // number of points of exact order n in Z/a x Z/b
  return &+[ MoebiusMu(n div e)*Gcd(a,e)*Gcd(b,e) : e in Divisors(n)];
end function;
cnt:=AssociativeArray();
for d in [1..D] do
  F:=GF(p^d); for n in Ns cat [-m: m in Ms] do cnt[<n,d>]:=0; end for;
  for j in F do
    E0:=EllipticCurveWithjInvariant(j);
    for T in Twists(E0) do
      N:=#T;
      rel:=[n: n in Ns | N mod n eq 0] cat [-m: m in Ms | N mod (2*m) eq 0];
      if #rel eq 0 then continue; end if;
      inv:=Invariants(AbelianGroup(T)); a:=#inv eq 2 select inv[1] else 1; b:=inv[#inv];
      aut:=(j eq 0 or j eq 1728) select #AutomorphismGroup(T) else 2;
      for n in rel do
        if n gt 0 then c:=exact(a,b,n);
        else m:=-n; // X_1(2,m): P of exact order m, Q in E[2] outside <P>; needs E[2] rational
          c:=(a mod 2 eq 0) select 2*exact(a,b,m) else 0; end if;
        assert c mod aut eq 0; cnt[<n,d>]+:=c div aut;
      end for;
    end for;
  end for;
end for;
for n in Ns cat [-m: m in Ms] do
  pts:=[cnt[<n,d>]: d in [1..D]];
  places:=[ (&+[MoebiusMu(d div e)*pts[e] : e in Divisors(d)]) div d : d in [1..D]];
  printf "PLACES %o p=%o: %o\n", n gt 0 select "X1("*IntegerToString(n)*")" else "X1(2,"*IntegerToString(-n)*")", p, places;
end for;
print "PLACECOUNTS_DONE"; quit;
