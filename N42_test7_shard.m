SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
if not assigned Shard then Shard:=0;elif Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;
if not assigned NShard then NShard:=4;elif Type(NShard) eq MonStgElt then NShard:=StringToInteger(NShard);end if;
if not assigned Limit then Limit:=0;elif Type(Limit) eq MonStgElt then Limit:=StringToInteger(Limit);end if;
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
X:=MDX1(42,GF(5));assert Genus(X) eq 25;
load "data/N42_degree7_representatives.m";
assert #N42D7Representatives eq 1860 and 0 le Shard and Shard lt NShard;
indices:=[i:i in [1..1860]|i mod NShard eq Shard];
print "N42_D7_SHARD_SETUP",Shard,NShard,#indices,Cputime(t0);
tested:=0;passing:=0;
for i in indices do
 if Limit gt 0 and tested ge Limit then break;end if;
 dd,def,coords:=Explode(N42D7Representatives[i]);assert dd eq 7;
 K:=ext<GF(5)|Polynomial(GF(5),def)>;
 ps:=Places(Curve(X)(K)![K!cs:cs in coords]);assert #ps eq 1 and Degree(ps[1]) eq 7;
 x:=ps[1];killed:=false;
 for q in [11,13,17] do
  E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
  y:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=q*P>);
  assert Degree(y) eq 7;
  T:=HeckeOperatorCommonZeros(X,q,x);
  A:=T-q*Divisor(y)-Divisor(x);
  if not IsPrincipal(A) then
   printf "N42_D7_SHARD_CASE index=%o q=%o nonprincipal=true time=%o\n",i,q,Cputime(t0);
   killed:=true;break;
  end if;
 end for;
 tested+:=1;
 if not killed then passing+:=1;print "N42_D7_SHARD_SURVIVOR",i;end if;
end for;
printf "N42_D7_SHARD_DONE shard=%o/%o tested=%o total=%o passing=%o complete=%o time=%o\n",Shard,NShard,tested,#indices,passing,tested eq #indices,Cputime(t0);
quit;
