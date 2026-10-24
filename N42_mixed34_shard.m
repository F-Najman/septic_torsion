/* New degree7 branches not supplied by the intrinsic lower-place tests:
   degree4 with larger cusp complement (q17 required), and mixed3+4. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(6*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("fast_diamond_hash_m2.m");
Attach("fast_diamond_hash_cyclic_bc.m");
X:=MDX1(42,GF(5));assert Genus(X) eq 25;
if Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;
assert Shard in [0..3];
load "data/N42_mixed34_fixed.m";
function LoadPlaces(data,X)
 out:=[];
 for row in data do
  d,def,co:=Explode(row);K:=ext<GF(5)|Polynomial(GF(5),def)>;
  ps:=Places(Curve(X)(K)![K!c:c in co]);assert #ps eq 1 and Degree(ps[1]) eq d;
  Append(~out,ps[1]);
 end for;
 assert #Seqset(out) eq #out;return out;
end function;
rs3:=LoadPlaces(N42Mixed34Fixed[1],X);ps4:=LoadPlaces(N42Mixed34Fixed[2],X);
assert #rs3 eq 2 and #ps4 eq 183;
print "N42_D7_MIXED_SHARD_SETUP",Shard,Cputime(t0);
function Aq(X,q,x)
 E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
 y:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=q*P>);
 return HeckeOperatorCommonZeros(X,q,x)-q*Divisor(y)-Divisor(x);
end function;
cache:=AssociativeArray();bad34:=[];idx:=0;
for i in [1..2] do
 for j in [1..183] do
  if (j-1) mod 4 ne Shard then continue;end if;
  idx+:=1;killed:=false;
  for q in [11,13,17] do
   for x in [rs3[i],ps4[j]] do
    key:=<q,x>;
    if not IsDefined(cache,key) then cache[key]:=Aq(X,q,x);end if;
   end for;
   if not IsPrincipal(cache[<q,rs3[i]>]+cache[<q,ps4[j]>]) then killed:=true;break;end if;
  end for;
  if not killed then Append(~bad34,<i,j>);end if;
  print "N42_D7_MIXED_SHARD_CASE",Shard,i,j,killed,Cputime(t0);
 end for;
end for;
expected:=2*#[j:j in [1..183]|(j-1) mod 4 eq Shard];assert idx eq expected;
printf "N42_D7_MIXED_SHARD_DONE shard=%o/4 tested=%o passing=%o complete=true time=%o\n",Shard,idx,#bad34,Cputime(t0);
print bad34;
quit;
