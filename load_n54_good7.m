/* Attach prime_degree_place.m before loading this function. */
function N54LoadGood7(X,filename)
 data:=eval Read(filename);def,records:=Explode(data);
 assert #records eq 828;
 S:=PolynomialRing(GF(5));f:=S!def;assert IsIrreducible(f) and Degree(f) eq 7;
 K<aa>:=ext<GF(5)|f>;FF:=FunctionField(Curve(X));out:=[];weighted:=0;
 for record in records do
  xycoeff,n,r,jcoeff:=Explode(record);
  assert #xycoeff eq 2 and n mod 54 eq 0 and r in [1,3];
  xy:=[Evaluate(S!co,aa):co in xycoeff];jj:=Evaluate(S!jcoeff,aa);
  assert Degree(MinimalPolynomial(jj,GF(5))) eq 7;
  // The inverse asserts the unique normalized branch, residue degree7,
  // and positive valuation of both exact coordinate-ideal generators.
  pp:=MDPlaceCommonZerosPrimeDegree([FF.1,FF.2],xy);
  assert Degree(pp) eq 7;Append(~out,pp);weighted+:=9*r;
 end for;
 assert #out eq 828 and #Seqset(out) eq 828 and weighted eq 11412;
 print "N54_GOOD7_LOADER_DONE representatives=828 weighted=11412 complete=true";
 return out;
end function;
