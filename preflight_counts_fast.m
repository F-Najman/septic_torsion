SetColumns(0);SetMemoryLimit(2*10^9);t0:=Cputime();
AttachSpec("gl2/magma.spec");
for tup in [<22,3>,<24,5>,<26,3>,<28,3>,<30,7>,<36,5>] do
 n,p:=Explode(tup);GG:=GL(2,Integers(n));
 HH:=sub<GG|[GG![1,0,2,1]] cat [GG![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1 and IsOdd(a)]>;
 zs:=GL2PointCounts(HH,p,7);
 ys:=[zs[d]-GL2RationalCuspCount(HH,p^d):d in [1..7]];
 assert &+GL2jCounts(HH,p) eq ys[1];closed:=[];
 for d in [1..7] do
  z:=ys[d]-&+([0] cat [e*closed[e]:e in [1..d-1]|d mod e eq 0]);
  assert z mod d eq 0;Append(~closed,z div d);
 end for;
 printf "DEG7_M2_PREFLIGHT N=%o p=%o genus=%o rational=%o closed=%o time=%o\n",n,p,GL2Genus(HH),ys,closed,Cputime(t0);
end for;
print "DEG7_M2_PREFLIGHT_DONE complete=true",Cputime(t0);
quit;
