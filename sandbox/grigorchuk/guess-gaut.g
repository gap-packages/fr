if true then
n := 7;
G := MakeGGroup(n);
A := x->Combine(2^(n-1),[(),Annihilate(x,2)]);
B := x->Combine(2^(n-1),[Annihilate(x,2),Annihilate(x,2)]);
a := G.1; b := G.2; c := G.3; d := G.4;
x := Comm(a,b); y := Comm(a,d);
N := ClosureGroup(G,[A(y),B(A(y)),B(B(A(y))),B(B(B(A(y))))]);
fi;
