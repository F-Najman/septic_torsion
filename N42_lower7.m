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
ps3:=NoncuspidalPlaces(X,3);ps4:=NoncuspidalPlaces(X,4);
assert #ps3 eq 4 and #ps4 eq 183;
rs3:=PlacesUpToDiamondCyclicHashBC(X,ps3);
rs4:=PlacesUpToDiamondCyclicHashBC(X,ps4);
assert #rs3 eq 2 and #rs4 eq 33;
print "N42_D7_LOWER_SETUP",Cputime(t0);
function Aq(X,q,x)
 E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
 y:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=q*P>);
 return HeckeOperatorCommonZeros(X,q,x)-q*Divisor(y)-Divisor(x);
end function;
bad4:=[];
for i in [1..33] do
 principal:=IsPrincipal(Aq(X,17,rs4[i]));
 print "N42_D7_GRADE4",i,principal,Cputime(t0);
 if principal then Append(~bad4,i);end if;
end for;
cache:=AssociativeArray();bad34:=[];idx:=0;
for i in [1..2] do
 for j in [1..183] do
  idx+:=1;killed:=false;
  for q in [11,13,17] do
   for x in [rs3[i],ps4[j]] do
    key:=<q,x>;
    if not IsDefined(cache,key) then cache[key]:=Aq(X,q,x);end if;
   end for;
   if not IsPrincipal(cache[<q,rs3[i]>]+cache[<q,ps4[j]>]) then killed:=true;break;end if;
  end for;
  if not killed then Append(~bad34,<i,j>);end if;
  print "N42_D7_MIXED34",i,j,killed,Cputime(t0);
 end for;
end for;
assert idx eq 366;
printf "N42_D7_LOWER_DONE grade4=33 survivors4=%o mixed34=366 survivors34=%o complete=true time=%o\n",#bad4,#bad34,Cputime(t0);
print bad4,bad34;
quit;
