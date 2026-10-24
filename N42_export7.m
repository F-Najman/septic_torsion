/* Export all degree7 representatives with exact field/place round trips.
   The canonical Tate hash checks coverage of every input closed place. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_diamond_hash_m2.m");
Attach("fast_diamond_hash_cyclic_bc.m");
X:=MDX1(42,GF(5));assert Genus(X) eq 25;
print "N42_D7_EXPORT_SETUP",Cputime(t0);
ps:=NoncuspidalPlaces(X,7);assert #ps eq 11160;
print "N42_D7_EXPORT_PLACES",#ps,Cputime(t0);
rs:=PlacesUpToDiamondCyclicHashBC(X,ps);
print "N42_D7_EXPORT_ORBITS",#rs,Cputime(t0);
data:=[* *];K:=GF(5^7);
for i in [1..#rs] do
 x:=rs[i];pt:=RepresentativePoint(x);vv:=MDM2ToCommonField(Eltseq(pt),K);
 def:=[Integers()!a:a in Eltseq(DefiningPolynomial(K,GF(5)))];
 coords:=[[Integers()!a:a in Eltseq(v,GF(5))]:v in vv];
 L:=ext<GF(5)|Polynomial(GF(5),def)>;
 xx:=Places(Curve(X)(L)![L!c:c in coords]);
 assert #xx eq 1 and xx[1] eq x and Degree(x) eq 7;
 Append(~data,[*7,def,coords*]);
 if i mod 100 eq 0 then print "N42_D7_EXPORT_PROGRESS",i,#rs,Cputime(t0);end if;
end for;
Write("data/N42_degree7_representatives.m",
 "N42D7Representatives:=" cat Sprint(data) cat ";" : Overwrite:=true);
printf "N42_D7_EXPORT_DONE places=11160 representatives=%o roundtrip=true complete=true time=%o\n",#rs,Cputime(t0);
quit;
