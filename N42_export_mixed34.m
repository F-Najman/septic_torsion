/* Export all degree7 representatives with exact field/place round trips.
   The canonical Tate hash checks coverage of every input closed place. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_diamond_hash_m2.m");
Attach("fast_diamond_hash_cyclic_bc.m");
X:=MDX1(42,GF(5));assert Genus(X) eq 25;
print "N42_D7_EXPORT_SETUP",Cputime(t0);
ps3:=NoncuspidalPlaces(X,3);ps4:=NoncuspidalPlaces(X,4);
assert #ps3 eq 4 and #ps4 eq 183;
rs3:=PlacesUpToDiamondCyclicHashBC(X,ps3);assert #rs3 eq 2;
records:=[* *];
for ps in [rs3,ps4] do
 data:=[* *];d:=Degree(ps[1]);K:=GF(5^d);
 for x in ps do
  pt:=RepresentativePoint(x);vv:=MDM2ToCommonField(Eltseq(pt),K);
  def:=[Integers()!a:a in Eltseq(DefiningPolynomial(K,GF(5)))];
  coords:=[[Integers()!a:a in Eltseq(v,GF(5))]:v in vv];
  L:=ext<GF(5)|Polynomial(GF(5),def)>;
  xx:=Places(Curve(X)(L)![L!c:c in coords]);
  assert #xx eq 1 and xx[1] eq x and Degree(x) eq d;
  Append(~data,[*d,def,coords*]);
 end for;
 Append(~records,data);
end for;
Write("data/N42_mixed34_fixed.m","N42Mixed34Fixed:=" cat Sprint(records) cat ";" : Overwrite:=true);
printf "N42_D7_MIXED_EXPORT_DONE cubic_reps=2 quartic_places=183 roundtrip=true complete=true time=%o\n",Cputime(t0);
quit;
