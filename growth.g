LoadPackage("automata");
growth := function(A,n)
    local v, i, j, m, s, r;
    A := UsefulAutomaton(A);
    s := A!.states;
    v := List([1..s],i->0);
    for i in InitialStatesOfAutomaton(A) do v[i] := 1; od;
    m := List([1..s],i->[]);
    r := UnderlyingMultiGraphOfAutomaton(A);
    for i in [1..s] do
        for j in r[i] do
            Add(m[j],i);
        od;
    od;
    r := [];
    for i in [1..n] do
        Add(r,Sum(v{FinalStatesOfAutomaton(A)}));
        v := List([1..s],i->Sum(v{m[i]}));
    od;
    return r;
end;

sphere := function(A,n,keep)
    local v, i, j, k, l, s, r, f, g;
    A := UsefulAutomaton(A);
    s := A!.states;
    v := List([1..s],i->[]);
    if IsList(AlphabetOfAutomaton(A)) then
        f := FreeMonoid(List(AlphabetOfAutomaton(A),x->[x]));
    else
        f := FreeMonoid(AlphabetOfAutomaton(A));
    fi;
    g := GeneratorsOfMonoid(f);
    for i in InitialStatesOfAutomaton(A) do v[i] := [One(f)]; od;
    for i in [1..n] do
        r := v;
        v := List([1..s],i->[]);
        for j in [1..A!.alphabet] do
            for k in [1..s] do
                for l in r[k] do
                    Add(v[A!.transitions[j][k]],l*g[j]);
                od;
            od;
        od;
    od;
    if keep then
        return Concatenation(v{FinalStatesOfAutomaton(A)});
    else
        return Concatenation(v{DifferenceLists([1..s],FinalStatesOfAutomaton(A))});
    fi;
end;

red := function(n)
    local t, g, d, i, N;
    t := [];
    N := 3^n;
    for g in [FabrykowskiGuptaGroup.1^-1,FabrykowskiGuptaGroup.2^-1,
            FabrykowskiGuptaGroup.1,FabrykowskiGuptaGroup.2] do
        d := Decomposition(g,n);
        Add(t,[]);
        for i in [1..N] do
            if IsOne(d[1][i]) then
                t[Length(t)][i] := i^d[2];
            fi;
        od;
    od;
    t[1][N+1] := [N+1,N+2]; t[1][N+2] := N+3; t[1][N+3] := N+3;
    t[2][N+1] := [N+1,N];   t[2][N] := N+3;   t[2][N+3] := N+3;
    t[3][N+1] := [N+1,N+2]; t[3][N+2] := N+3; t[3][N+3] := N+3;
    t[4][N+1] := [N+1,N];   t[4][N] := N+3;   t[4][N+3] := N+3;
    return Automaton("nondet",N+3,"ATat",t,[N+1],[N+3]);
end;

Letter := RationalExpression("aUtUAUT");
Cancellation := RationalExpression("((tUT)(tUT))U((aUA)(aUA))");
Word := StarRatExp(Letter);
LStar := RatExpToAut(Word);
Tt := Automaton("det",2,"ATat",[[],[2],[],[2]],[1],[2]);
#Tt := RationalExpression("tUT","ATat");
NonReduced := ProductRatExp(ProductRatExp(Word,Cancellation),Word);
NR := RatExpToAut(NonReduced);
#for i in [1..3] do
#    NR := UnionAutomata(NR,red(i));
#od;
m := FRMachine(["A","T","a","t"],[[[],[],[]],[[1],[],[2]],[[],[],[]],[[3],[],[4]]],[(1,3,2),(),(1,2,3),()]);
A := FRElement(m,1); T := FRElement(m,2);
a := FRElement(m,3); t := FRElement(m,4);
