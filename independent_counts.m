SetColumns(0);
SetMemoryLimit(10^9);
AttachSpec("gl2/magma.spec");
for t in [<1,32,3>,<1,34,3>,<1,35,3>,<1,38,3>,<1,40,3>,<1,50,3>,<1,39,5>] do
 m,n,p:=Explode(t);
 G:=GL(2,Integers(n));
 H:=sub<G | [G![1,0,m,1]] cat
  [G![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1 and (a-1) mod m eq 0]>;
 ycnt:=[&+GL2jCounts(H,p^d):d in [1..7]];
 a:=[];
 for d in [1..7] do
  z:=ycnt[d]-&+([0] cat [e*a[e]:e in [1..d-1]|d mod e eq 0]);
  assert z mod d eq 0;
  Append(~a,z div d);
 end for;
 printf "DEGREE7_INDEPENDENT_COUNT m=%o n=%o p=%o genus=%o rational_point_counts=%o closed_place_counts=%o\n",
  m,n,p,GL2Genus(H),ycnt,a;
end for;
print "DEGREE7_INDEPENDENT_COUNT_DONE";
quit;
