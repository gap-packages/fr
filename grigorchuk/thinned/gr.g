RequirePackage("treegp");

n := 7;
p := 2;
K := GF(p);

G := Image(EpimorphismPermGroupFrGroup(GrigorchukGroup,n));
S := List(GeneratorsOfGroup(G),g->PermutationMat(g,2^n,K));

oldgr := function(l)
  local V, W;

  V := [VectorSpace(K,[S[1]^2])];
  while Size(V) < l do
#    W := V[Size(V)];
#    for p in Cartesian(Basis(W),S) do
#      W := ClosureLeftModule(W,p[1]*p[2]);
#    od;
    W := Basis(V[Size(V)]);
    Add(V,VectorSpace(K,Concatenation(W,List(Cartesian(W,S),p->p[1]*p[2]))));
  od;
  return V;
end;

gr := function(l)
  local B, V, b, q;

  b := [S[1]^2];
  B := [b];
  while Size(B) < l do
    V := VectorSpace(K,Concatenation(b,List(Cartesian(B[Size(B)],S),p->p[1]*p[2])));
    q := Basis(V);
    Add(B,q{[Size(b)+1..Size(q)]});
    b := q;
    Print("---> ",Size(b),"\n");
  od;
  return B;
end;
