/* Exact coordinate ideals in prime-degree finite fields, with no new
   subfield objects or generic field embeddings. */
intrinsic MDTwoGeneratorsPrimeDegree(x::FldFinElt,y::FldFinElt) -> RngMPolElt,RngMPolElt
{ Generators of the closed coordinate point, when the pair has prime degree }
 K:=Parent(x);require Parent(y) cmpeq K : "Require a common explicit field parent";
 F:=PrimeField(K);d:=Degree(K,F);require IsPrime(d) : "Require prime degree";
 fx:=MinimalPolynomial(x,F);fy:=MinimalPolynomial(y,F);
 assert Degree(fx) in [1,d] and Degree(fy) in [1,d];
 require Degree(fx) eq d or Degree(fy) eq d : "The pair must generate the field";
 R<u,v>:=PolynomialRing(F,2);
 if Degree(fx) eq d then
  M:=Matrix(F,[Eltseq(x^i,F):i in [0..d-1]]);assert Rank(M) eq d;
  coeff:=Solution(M,Vector(F,Eltseq(y,F)));
  assert &+[K!coeff[i+1]*x^i:i in [0..d-1]] eq y;
  f:=Evaluate(fx,u);g:=v-&+[coeff[i+1]*u^i:i in [0..d-1]];
 else
  M:=Matrix(F,[Eltseq(y^i,F):i in [0..d-1]]);assert Rank(M) eq d;
  coeff:=Solution(M,Vector(F,Eltseq(x,F)));
  assert &+[K!coeff[i+1]*y^i:i in [0..d-1]] eq x;
  f:=Evaluate(fy,v);g:=u-&+[coeff[i+1]*v^i:i in [0..d-1]];
 end if;
 assert Evaluate(f,[x,y]) eq 0 and Evaluate(g,[x,y]) eq 0;
 return f,g;
end intrinsic;

intrinsic MDPlaceCommonZerosPrimeDegree(coordinates::[FldFunFracSchElt],values::[FldFinElt]) -> PlcCrvElt
{ Recover the normalized place using an exact prime-degree coordinate ideal }
 require #coordinates eq 2 and #values eq 2 : "Require two coordinates and values";
 f,g:=MDTwoGeneratorsPrimeDegree(values[1],values[2]);
 f1:=Evaluate(f,coordinates);g1:=Evaluate(g,coordinates);
 FF:=Parent(coordinates[1]);C:=Curve(FF);AF,ma:=AlgorithmicFunctionField(FF);
 O:=MaximalOrderFinite(AF);I:=ideal<O|ma(f1),ma(g1)>;
 fs:=Factorization(I);ys:=[CurvePlace(C,Place(z[1])):z in fs|z[2] gt 0];
 if #ys eq 0 then
  Oi:=MaximalOrderInfinite(AF);Ii:=ideal<Oi|ma(f1),ma(g1)>;
  fsi:=Factorization(Ii);ys:=[CurvePlace(C,Place(z[1])):z in fsi|z[2] gt 0];
 end if;
 assert #ys eq 1;y:=ys[1];
 assert Degree(y) eq Degree(Parent(values[1]),PrimeField(Parent(values[1])));
 assert Valuation(f1,y) gt 0 and Valuation(g1,y) gt 0;
 return y;
end intrinsic;

intrinsic ModuliPointCyclicPrimeDegree(X::MDCrvMod1,E::CrvEll,Ls::Rec) -> PlcCrvElt
{ Cyclic inverse for a full pair of prime residue degree, avoiding subfields }
 xy:=TateNormalForm_xy(E,Ls`P);FF:=FunctionField(Curve(X));
 return MDPlaceCommonZerosPrimeDegree([FF.1,FF.2],xy);
end intrinsic;
