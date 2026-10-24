/* Explicitly restrict Decomposition candidates to the selected coefficient
   place. Magma2.30-1's infinite-place path returns all places above the
   underlying rational infinity, even when an intermediate place is given.
   See probe_relative_infinite_decomposition.m/.log for a 0.050CPU example. */
intrinsic DecompositionAboveSelectedPlace(L::FldFun,P::PlcFunElt : Unramified:=false) -> SeqEnum
{All places above the exact immediate-coefficient-field place P, with a degree audit}
 K:=FunctionField(P);
 require BaseField(L) cmpeq K : "P must belong to the exact immediate coefficient field";
 f,g:=TwoGenerators(P);
 assert Valuation(f,P) gt 0 and Valuation(g,P) gt 0;
 candidates:=Decomposition(L,P);
 above:=[Q:Q in candidates|Valuation(L!f,Q) gt 0 and Valuation(L!g,Q) gt 0];
 assert #above gt 0;
 u:=UniformizingElement(P);assert Valuation(u,P) eq 1;
 es:=[Valuation(L!u,Q):Q in above];assert &and[e gt 0:e in es];
 assert &+[es[i]*Degree(above[i]):i in [1..#above]] eq Degree(L)*Degree(P);
 if Unramified then assert &and[e eq 1:e in es];end if;
 return above;
end intrinsic;
