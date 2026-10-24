/* Fresh complete pure7 input, independent of every old degree7 index. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
if Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;
assert Shard in [0,1];
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("prime_degree_place.m");
X:=MDX11(2,28,GF(3));cs:=Cusps(X);
def,records:=Explode(eval Read("data/N28_D7_hecke_input.txt"));
assert #records eq 52;R:=PolynomialRing(GF(3));f:=R!def;
assert IsIrreducible(f) and Degree(f) eq 7;K<aa>:=ext<GF(3)|f>;
tested:=0;survivors:=[];
for i in [1..52] do
 if (i-1) mod 2 ne Shard then continue;end if;
 bc:=[Evaluate(R!co,aa):co in records[i]];
 x:=MDPlaceCommonZerosPrimeDegree(X`_coordinates,bc);assert Degree(x) eq 7;
 killed:=false;
 for q in [5,11] do
  A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));
  assert Degree(A) eq 0;principal:=IsPrincipal(A);
  print "N28_D7_FRESH_HECKE_CASE",Shard,i,q,principal,Cputime(t0);
  if not principal then killed:=true;break;end if;
 end for;
 tested+:=1;
 if not killed then Append(~survivors,i);end if;
end for;
assert tested eq 26;
printf "N28_D7_FRESH_HECKE_DONE shard=%o tested=%o remaining=%o complete=true time=%o\n",Shard,tested,#survivors,Cputime(t0);
print "N28_D7_FRESH_HECKE_SURVIVORS",survivors;quit;
