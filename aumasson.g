RndTreeAut := function(br)
    local trans, out, done, todo, i, j, t;
    trans := [];
    out := [];
    todo := 1; done := 0;
    while todo>done do
        Add(out,Random([(),(1,2)]));
        t := [];
        for j in [1,2] do
            if Random([1..DenominatorRat(br)])<=NumeratorRat(br) then
                todo := todo+1;
                Add(t,todo);
            else
                Add(t,1);
            fi;
        od;
        Add(trans,t);
        done := done+1;
    od;
    return MealyElement(trans,out,1);
end;

RndTreeAut := function(n)
    local e, i, j;
    e := FR_LOCAL.RANDOMFINITARY(FullBinaryGroup,n);
    for i in e!.transitions do for j in [1..2] do
        if i[j]=Length(e!.transitions) then i[j] := 1; fi;
    od; od;
    return Minimized(e);
end;
