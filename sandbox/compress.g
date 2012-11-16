inc := function(G,n)
    local S, s, g, h, gens, todo;
    gens := List(Difference(Nucleus(G),[One(G)]),g->[g,Length(Germs(g))]);
    S := List([1..n],i->[]);
    todo := [[One(G),0]];
    while not IsEmpty(todo) do
        g := Remove(todo,1);
        if g[1] in S[1+g[2]] then continue; fi;
        AddSet(S[1+g[2]],g[1]);
        for s in gens do
            if g[2]+s[2]<n then
                h := g[1]*s[1];
                if Length(Germs(h))=g[2]+s[2] then
                    Add(todo,[h,g[2]+s[2]]);
                fi;
            fi;
        od;
    od;
    return S;
end;

# grig. gp: 836 elts
# [ 2, 12, 34, 72, 142, 144, 144, 144, 142 ]

# 0|1: [ 2, 4, 4, ... ]
# 00|1: [ 8, 64, 320, 1280, 4288, 12288, ? ]
# 0|11: [ 2, 12, 34, 72, 142, 144, 144, ... ]
# 0|111: [ 2, 28, 178, 784, 2976, 6272, 12544, ? ]
# 0|01: [ 2, 8, 22, 52, 110, 228, 430, 740, 1144, 1652, 2304, 3184, 4382, 5980, 8008, 10388, 12926, ??? ]
# 0|011: [ 2, 12, 44, 128, 326, 760, 1650, ? ]
# 0|0101: [ 2, 16, 70, 236, 676, 1708, ? ]
