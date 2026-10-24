/* Independent Newton-sum check from the modular-symbol Frobenius
   polynomial against the complete function-field place enumeration. */
assert Degree(fp) eq 42;trace33:=[];closed33:=[];
for n33 in [1..7] do
 sn33:=-n33*Coefficient(fp,42-n33);
 for j33 in [1..n33-1] do sn33-:=Coefficient(fp,42-j33)*trace33[n33-j33];end for;
 Append(~trace33,Integers()!sn33);count33:=2^n33+1-Integers()!sn33;
 for d33 in Divisors(n33) do
  if d33 lt n33 then count33-:=d33*closed33[d33];end if;
 end for;
 assert count33 mod n33 eq 0;Append(~closed33,count33 div n33);
 assert closed33[n33] eq #places[n33]+#[c:c in cs|Degree(c) eq n33];
end for;
print "C33_FROBENIUS_COUNTS_VERIFIED noncuspidal",[#ss:ss in places],"all_places",closed33;
