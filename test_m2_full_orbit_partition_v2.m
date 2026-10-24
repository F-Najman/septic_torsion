/* Independent carried-point inverse versus direct (b,c) orbit hashing. */
SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
if not assigned N then N:=26;elif Type(N) eq MonStgElt then N:=StringToInteger(N);end if;
if not assigned d then d:=5;elif Type(d) eq MonStgElt then d:=StringToInteger(d);end if;
assert N in [26,28] and d in [5,6,7];
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("fast_diamond_hash_m2.m");
X:=MDX11(2,N,GF(3));S:=NoncuspidalPlaces(X,d);
rs:=PlacesUpToDiamondM2HashDirect(X,S);
expected:=N eq 26 select [6,26,48] else [8,16,52];
assert #rs eq expected[d-4];
print "M2_ORBIT_REFERENCE_SETUP",N,d,#S,#rs,Cputime(t0);
K:=Parent(Evaluate(X`_coordinates[1],S[1]));
LR:=recformat<P:PtEll,Q:PtEll>;
covered:={S[1]};Exclude(~covered,S[1]);all:=Seqset(S);
for i in [1..#rs] do
 bc:=MDM2CanonicalPair(X,rs[i],K);b,c:=Explode(bc);
 E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P2:=E![0,0];Q:=E![b,b];
 orbit:={rs[i]};
 for a in [1..N div 2] do
  if GCD(a,N) ne 1 then continue;end if;
  L:=rec<LR|P:=P2,Q:=a*Q>;
  Include(~orbit,ModuliPointFast(X,E,L));
 end for;
 print "M2_ORBIT_REFERENCE_SINGLE",N,d,i,#orbit,Cputime(t0);
 assert orbit subset all and IsEmpty(orbit meet covered);
 covered join:=orbit;
 if i mod 10 eq 0 then print "M2_ORBIT_REFERENCE_PROGRESS",N,d,i,#rs,#covered,Cputime(t0);end if;
end for;
assert covered eq all;
printf "M2_ORBIT_REFERENCE_DONE N=%o degree=%o places=%o orbits=%o complete=true time=%o\n",N,d,#S,#rs,Cputime(t0);
quit;
