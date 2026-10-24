/* Tiny exact regression: relative infinite-place decomposition in a tower.
   K:a^2+a=1/x; L:b^2+b=a. Both fields are rational, and x=infinity
   has a=0 (split in L) and a=1 (inert in L). */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(5*10^8);t0:=Cputime();
F<x>:=FunctionField(GF(2));R<u>:=PolynomialRing(F);
K<a>:=FunctionField(u^2+u+1/x);S<v>:=PolynomialRing(K);
L<b>:=FunctionField(v^2+v+a);
assert Genus(K) eq 0 and Genus(L) eq 0;
cs:=Poles(K!x);assert #cs eq 2 and &and[Degree(P) eq 1:P in cs];
for P in cs do
 aa:=Evaluate(a,P);assert aa in {GF(2)!0,GF(2)!1};
 raw:=Decomposition(L,P);
 f,g:=TwoGenerators(P);
 filtered:=[Q:Q in raw|Valuation(L!f,Q) gt 0 and Valuation(L!g,Q) gt 0];
 printf "RELATIVE_INFINITY a=%o raw_degrees=%o filtered_degrees=%o\n",aa,[Degree(Q):Q in raw],[Degree(Q):Q in filtered];
 expected:=aa eq 0 select [1,1] else [2];
 assert Sort([Degree(Q):Q in filtered]) eq expected;
 assert &+[Degree(Q):Q in filtered] eq 2;
end for;
print "RELATIVE_INFINITY_REGRESSION_DONE complete=true",Cputime(t0);
quit;
