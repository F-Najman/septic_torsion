/* Independent Frobenius/Newton count, separate from place enumeration. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);t0:=Cputime();
A:=JH(56,[29]);assert Dimension(A) eq 25;
f:=FrobeniusPolynomial(A,3);assert Degree(f) eq 50;
co:=[Integers()!Coefficient(f,i):i in [0..50]];
assert co[51] eq 1;
ss:=[];closed:=[];
for r in [1..7] do
 s:=-r*co[51-r];
 for i in [1..r-1] do s-:=co[51-i]*ss[r-i];end for;
 Append(~ss,s);
 nr:=3^r+1-s-&+[Integers()|d*closed[d]:d in Divisors(r)|d lt r];
 assert nr mod r eq 0;Append(~closed,nr div r);
end for;
cusps:=[12,6,4,0,0,2,0];nc:=[closed[i]-cusps[i]:i in [1..7]];
assert nc eq [0,0,0,0,48,64,312];
Write("data/N28_D7_frobenius_counts.txt",Sprint([*co,ss,closed,cusps,nc*]) : Overwrite:=true);
print "N28_D7_FROBENIUS_COUNTS",nc,"complete=true",Cputime(t0);quit;
