/* Exact kernel-orbit Hecke-2 on an odd-level curve in odd characteristic. */
intrinsic HeckeTwoCommonZeros(X::MDCrvMod1,x::PlcCrvElt) -> DivCrvElt
{ Exact T_2 with common-zero place recovery }
    E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;K:=BaseRing(E);
    require IsOdd(Characteristic(K)) and IsOdd(Level(X)):
        "Require odd level and odd residue characteristic";
    b2,b4,b6,b8:=Explode(bInvariants(E));
    R<z>:=PolynomialRing(K);
    f:=4*z^3+b2*z^2+2*b4*z+b6;
    assert Degree(f) eq 3 and GCD(f,Derivative(f)) eq 1;
    answer:=DivisorGroup(Curve(X))!0;sumorbits:=0;
    for fac in Factorization(f) do
        h:=fac[1];e:=Degree(h);assert fac[2] eq 1;
        if e eq 1 then L:=K;aa:=-Coefficient(h,0)/Coefficient(h,1);
        else L<aa>:=ext<K|h>;end if;
        EL:=EllipticCurve([L!a:a in aInvariants(E)]);
        RL<w>:=PolynomialRing(L);
        EG,phi:=IsogenyFromKernel(EL,w-aa);
        PG:=phi(EL![L!P[1],L!P[2]]);
        y:=ModuliPointCyclicCommonZeros(X,EG,rec<recformat<P:PtEll>|P:=PG>);
        mult:=e*Degree(x)/Degree(y);assert Denominator(mult) eq 1;
        answer+:=(Integers()!mult)*Divisor(y);sumorbits+:=e;
    end for;
    assert sumorbits eq 3 and Degree(answer) eq 3*Degree(x);
    return answer;
end intrinsic;
