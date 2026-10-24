/* Actual presentations include the free degree generator as a final
   zero column; there is no redundant zero relation row. */
RG36:=RelationMatrix(G);RY36:=RelationMatrix(GY36);
ig36:=Invariants(G);iyg36:=Invariants(GY36);
assert Ngens(G) eq #ig36 and Ngens(GY36) eq #iyg36;
assert ig36[#ig36] eq 0 and iyg36[#iyg36] eq 0;
eg36:=HorizontalJoin(DiagonalMatrix(Integers(),Prune(ig36)),ZeroMatrix(Integers(),#ig36-1,1));
ey36:=HorizontalJoin(DiagonalMatrix(Integers(),Prune(iyg36)),ZeroMatrix(Integers(),#iyg36-1,1));
assert RG36 eq eg36 and RY36 eq ey36;
lg36:=Open("logs/x1_36/coordinate_presentations.log","w");
fprintf lg36,"C36_SOURCE_GENERATORS ngens=%o invariants=%o\n",Ngens(G),ig36;
fprintf lg36,"C36_SOURCE_ACTUAL_RELATIONS\n%o\n",RG36;
fprintf lg36,"C36_TARGET_GENERATORS ngens=%o invariants=%o\n",Ngens(GY36),iyg36;
fprintf lg36,"C36_TARGET_ACTUAL_RELATIONS\n%o\n",RY36;
fprintf lg36,"C36_COORDINATE_PRESENTATIONS_FINITE_DIAGONAL_VERIFIED complete=true\n";
Flush(lg36);delete lg36;
cd36:=Open("data/coordinate_presentations_36_data.py","w");
fprintf cd36,"source_relations = %o\n",[Eltseq(r):r in Rows(RG36)];
fprintf cd36,"target_relations = %o\n",[Eltseq(r):r in Rows(RY36)];
fprintf cd36,"source_ngens = %o\n",Ngens(G);
fprintf cd36,"target_ngens = %o\n",Ngens(GY36);
delete cd36;
print "C36_COORDINATE_PRESENTATIONS_FINITE_DIAGONAL_VERIFIED complete=true";
