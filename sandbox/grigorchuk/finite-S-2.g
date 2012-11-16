RequirePackage("treegp");

p := 3; K := GF(p); o := One(K); z := Zero(K);

if not IsBound(A) then
A := 0; B := 0; C := 0; D := 0; # shut up warnings
a := 0; b := 0; c := 0; d := 0;
fi;

FiniteThinnedAlgebra := function(g,F)
  local n;
  n := Maximum(List(GeneratorsOfGroup(g),LargestMovedPoint));
#  F := GF(PrimePGroup(g));
  return RingByGenerators(List(GeneratorsOfGroup(g),x->PermutationMat(x,n,F)));
end;

ThinnedAugmentationIdeal := function(a)
  return Ideal(a,List(GeneratorsOfAlgebra(a),x->x-One(a)));
end;

#
# product of left-ideal by right-ideal, giving a two-sided ideal
# (the only case where it is easy to compute)
#
ProductOfIdeals := function(i,j)
  local l, r, fixup;
  fixup := function(i,h_aroi,hXaroi,a_i,Xaroi)
    if h_aroi(i) then return i;
    elif hXaroi(i) then return a_i(Xaroi(i),i);
    else return a_i(i,i);
    fi;
  end;
  l := fixup(i,HasLeftActingRingOfIdeal,HasRightActingRingOfIdeal,AsLeftIdeal,RightActingRingOfIdeal);
  r := fixup(j,HasRightActingRingOfIdeal,HasLeftActingRingOfIdeal,AsRightIdeal,LeftActingRingOfIdeal);
  return Ideal(LeftActingRingOfIdeal(l),
    Union(List(GeneratorsOfLeftIdeal(l),x->List(GeneratorsOfRightIdeal(r),y->x*y))));
end;

ProdA := function(i,j)
  return Ideal(LeftActingRingOfIdeal(i),
    Union(List(GeneratorsOfIdeal(i),x->List(GeneratorsOfIdeal(j),y->x*y))));
end;

ProdAA := function(i,j)
  return Ideal(LeftActingRingOfIdeal(i),
    Union(List(GeneratorsOfAlgebra(i),x->List(GeneratorsOfAlgebra(j),y->x*y))));
end;

ProdIA := function(i,j)
  return Ideal(LeftActingRingOfIdeal(i),
    Union(List(GeneratorsOfIdeal(i),x->List(GeneratorsOfAlgebra(j),y->x*y))));
end;

LieProd := function(i,j)
  return Ideal(LeftActingRingOfIdeal(i),
    Union(List(GeneratorsOfIdeal(i),x->List(GeneratorsOfIdeal(j),y->x*y-y*x))));
end;

BranchingIdeal := function(A)
  local g;
  g := GeneratorsOfAlgebra(A);
  if p=2 then
    return Ideal(A,[g[2]*g[1]-g[1]*g[2],(One(A)+g[1]*g[2])*(One(A)+g[1]),
      (g[4]+g[1]*g[3])*(One(A)+g[1])]);
  else
    return Ideal(A,[g[2]*g[1]-g[1]*g[2]]);
  fi;
end;

DSA := function(A)
  local L, k, i;
  L := [AsTwoSidedIdeal(A,A)];
  i := ThinnedAugmentationIdeal(A); k := i;
  while not IsTrivial(k) do
    Add(L,k);
    k := ProdA(k,i);
  od;
  return L;
end;

LSA := function(A)
  local L, k, i;
  L := [AsTwoSidedIdeal(A,A)];
  i := ThinnedAugmentationIdeal(A); k := i;
  while not IsTrivial(k) do
    Add(L,k);
    k := LieProd(k,i);
  od;
  return L;
end;

_1000 := x->KroneckerProduct([[1,0],[0,0]]*o,x);
_0100 := x->KroneckerProduct([[0,1],[0,0]]*o,x);
_0010 := x->KroneckerProduct([[0,0],[1,0]]*o,x);
_0001 := x->KroneckerProduct([[0,0],[0,1]]*o,x);
_1100 := x->KroneckerProduct([[1,1],[0,0]]*o,x);
_0011 := x->KroneckerProduct([[0,0],[1,1]]*o,x);
_1001 := x->KroneckerProduct([[1,0],[0,1]]*o,x);
_0110 := x->KroneckerProduct([[0,1],[1,0]]*o,x);
_1111 := x->KroneckerProduct([[1,1],[1,1]]*o,x);

if false and not IsBound(A5) then
  A5 := FiniteThinnedAlgebra(Image(EpimorphismPermGroupFrGroup(GrigorchukGroup,5)));
  g5 := GeneratorsOfAlgebra(A5);
  I5 := ThinnedAugmentationIdeal(A5);
  A4 := FiniteThinnedAlgebra(Image(EpimorphismPermGroupFrGroup(GrigorchukGroup,4)));
  g4 := GeneratorsOfAlgebra(A4);
  I4 := ThinnedAugmentationIdeal(A4);
  A3 := FiniteThinnedAlgebra(Image(EpimorphismPermGroupFrGroup(GrigorchukGroup,3)));
  g3 := GeneratorsOfAlgebra(A3);
  I3 := ThinnedAugmentationIdeal(A3);
# M3 := MatrixAlgebra(A3,2);
  P := x->x{[1..Size(x)/2]}{[1..Size(x)/2]};
  M3 := Ring(Union([g4[1]],List(g3,_1000)));
  A3i := Ring(List(g3,_1000));
  M4 := Ring(Union([g5[1]],List(g4,_1000)));
  A4i := Ring(List(g4,_1000));

  D3 := DSA(A3);
  D4 := DSA(A4);
  D5 := DSA(A5);

  i1000 := x->Ideal(A5,List(GeneratorsOfIdeal(x),_1000));
  i1100 := x->Ideal(A5,List(GeneratorsOfIdeal(x),_1100));
  i1001 := x->Ideal(A5,List(GeneratorsOfIdeal(x),_1001));
  i1111 := x->Ideal(A5,List(GeneratorsOfIdeal(x),_1111));

  D3_1000 := List(D3,i1000);
  D3_1100 := List(D3,i1100);
  D3_1001 := List(D3,i1001);
  D3_1111 := List(D3,i1111);
  D4_1000 := List([1..Size(D4)],
    function(n) if n>=4 then return i1000(D4[n]); else return A5; fi; end);
  D4_1100 := List([1..Size(D4)],
    function(n) if n>=4 then return i1100(D4[n]); else return A5; fi; end);
  D4_1001 := List([1..Size(D4)],
    function(n) if n>=4 then return i1001(D4[n]); else return A5; fi; end);
  D4_1111 := List([1..Size(D4)],
    function(n) if n>=4 then return i1111(D4[n]); else return A5; fi; end);
fi;

Union2Spaces := function(x,y)
  return VectorSpace(GF(2),Union(GeneratorsOfVectorSpace(x),GeneratorsOfVectorSpace(y)));
end;

#
# codimensions:
#
#0 A
#1 D
#2 |
#3 |
#4 D^2
#5 |
#6 K
#7 |
#8 |
#9 D^3
#
# D^n/D^(n+1): [ 1,  3,  4,  5,  6,  8, 10, 11, 12, 14,
#               16, 18, 20, 21, 22, 23, 24, 26, 28, 30,
#               32, 34, 36, 38, 40, 41, ... ]
#
# A  /D   = 1
# D  /D^2 = [A,B,C] (since D=B+C)
# D^2/D^3 = [AB,AC,BA,CA] (since A^2=B^2=C^2=BC=CB=0)
# D^3/D^4 = [ABA,ACA,BAB,BAC,CAB] (since CAC+CAB+BAC+BAB=0)
# D^4/D^5 = [ABAB,ABAC,ACAB,BABA,BACA,CABA]
# D^5/D^6 = [ABABA,ABACA,ACABA,BABAB,BABAC,BACAB,CABAB,CABAC]
# D^6/D^7 = [ABABAB,ABABAC,ABACAB,ACABAB,ACABAC,
#            BABABA,BABACA,BACABA,CABABA,CABACA]
# D^7/D^8 = [ABABABA,ABABACA,ABACABA,ACABABA,ACABACA,
#            BABABAC,BABACAB,BACABAB,BACABAC,CABABAB,CABACAB]
# (since BACABAC+CABACAB+BABABAB=0
#    and BABABAC+BACABAC+CABABAB+CABABAC+CABACAB=0)
# D^8/D^9 = [ABABABAC,ABABACAB,ABACABAB,ABACABAC,ACABABAB,ACABACAB,
#            BABABACA,BABACABA,BACABABA,BACABACA,CABABABA,CABACABA]
# conjecture:
# dim D^n/D^(n+1) = 2n-2^(x-1) if       2^x<=n<=3.2^(x-1)
#                 =  n+2^x     if 3.2^(x-1)<=n<=2^(x+1)

gen := function(s,ins,out)
  return Filtered(List(Cartesian(s,[A,B,C]),x->x[1]*x[2]),x->x in ins and not x in out);
end;

words := function(n)
  if n=1 then return [b,c,d];
  else return List(Cartesian(words(n-1),[b,c,d]),w->w[1]*a*w[2]);
  fi;
end;
semi := function(n)
  local s;
  s := FreeSemigroup("A","B","C","D","Z");
  a := GeneratorsOfSemigroup(s)[1];
  b := GeneratorsOfSemigroup(s)[2];
  c := GeneratorsOfSemigroup(s)[3];
  d := GeneratorsOfSemigroup(s)[4];
  z := GeneratorsOfSemigroup(s)[5];
  s := s / Union([[z^2,z],[z*a,z],[z*b,z],[z*c,z],[z*d,z],
	[a*z,z],[b*z,z],[c*z,z],[d*z,z],
	[a^2,z],[b^2,z],[c^2,z],[d^2,z],
	[b*c,z],[b*d,z],[c*d,z],[c*b,z],[d*b,z],[d*c,z],
	[d*a*d,z],[c*a*c*a*c*a*c,z],
	[d*a*b*a*b*a*d,z],[d*a*b*a*c*a*d,z],[d*a*c*a*b*a*d,z],[d*a*c*a*c*a*d,z]
	],
	List(words(n),w->[w*a,z]),List(words(n),w->[a*w,z]));
  a := GeneratorsOfSemigroup(s)[1];
  b := GeneratorsOfSemigroup(s)[2];
  c := GeneratorsOfSemigroup(s)[3];
  d := GeneratorsOfSemigroup(s)[4];
  z := GeneratorsOfSemigroup(s)[5];
  return s;
end;
if false then
  a := FreeAssociativeAlgebraWithOne(GF(2),4);
  A := a.1; B := a.2; C := a.3; D := a.4;
  i := Ideal(a,[A^2,B^2,C^2,D^2,B*C,C*B,B*D,D*B,C*D,D*C,B+C+D,
	D*A*D,C*A*C*A*C*A*C,D*A*B*A*B*A*D]);
  a := a/i;
  A := GeneratorsOfAlgebra(a)[1];
  B := GeneratorsOfAlgebra(a)[2];
  C := GeneratorsOfAlgebra(a)[3];
  D := GeneratorsOfAlgebra(a)[4];
fi;
WORDS := function(n)
  if n=1 then return [B,C,D];
  else return List(Cartesian(WORDS(n-1),[B,C,D]),w->w[1]*A*w[2]);
  fi;
end;
PP := function(m)
  local f;
  f := function(z) if IsZero(z) then return "  "; else return "* "; fi; end;
  Print(Concatenation(List(m,x->Concatenation(Concatenation(List(x,f)),"\n"))));
end;
AlgIndex := function(g)
  return Size(Semigroup(g));
end;

ggroup := n->Image(EpimorphismPermGroupFrGroup(GrigorchukGroup,n));
gens := G->List(GeneratorsOfGroup(G),g->PermutationMat(g,LargestMovedPoint(G),K));
gengrowth := function(S,B0,l)
  local B, V, b, q;

  b := B0;
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
alggrowth := function(S,l)
  return gengrowth(S,[S[1]^2],l);
end;
idgrowth := function(S,l)
  local x;
  x := S[1]*S[2]-S[2]*S[1];
  return gengrowth(S,[x,S[4]*x,S[1]*S[4]*x],l);
end;

_xxxx := V->VectorSpace(K,List(Basis(V),_1111));
_xooo := V->VectorSpace(K,Concatenation(List(Basis(V),_1000),List(Basis(V),_0100),List(Basis(V),_0010),List(Basis(V),_0001)));
_xoox := V->VectorSpace(K,Concatenation(List(Basis(V),_1001),List(Basis(V),_0110)));
_xxoo := V->VectorSpace(K,Concatenation(List(Basis(V),_1100),List(Basis(V),_0011)));

ns := 9;
nl := 9;

if false then
Ss := gens(ggroup(7));
Bs := alggrowth(Ss,ns);
Vs := VectorSpace(K,Concatenation(Bs));
Ws := VectorSpace(K,Concatenation(idgrowth(Ss,8)));
Fs := List([1..ns],i->VectorSpace(K,Concatenation(Bs{[1..i]})));
Es := List(Fs,v->Intersection(v,Ws));

Sl := gens(ggroup(8));
Bl := alggrowth(Sl,nl);
Vl := VectorSpace(K,Concatenation(Bl));
Wl := VectorSpace(K,Concatenation(idgrowth(Sl,8)));
Xl := Intersection(Wl,_xooo(Ws));
Fl := List([1..nl],i->VectorSpace(K,Concatenation(Bl{[1..i]})));
El := List(Fl,v->Intersection(v,Wl));
Dl := List(Fl,v->Intersection(v,Xl));
fi;
f := function(n)
  local k;

  if n<=6 then k := [4,6,8,10,13,16]; return k[n]; fi;

  k := 2^LogInt(n,2);
  if   k<=    n and n <= 5/4*k then return 4*n-3/2*k;
  elif k<=5/4*n and n <= 3/2*k then return 3*n-1/4*k;
  elif k<=3/2*n and n <= 7/4*k then return  n+11/4*k;
  elif k<=7/4*n and n <=   2*k then return 2*n+k;
  fi;
end;

K := GF(5);
S := gens(ggroup(7));
V := VectorSpace(K,Concatenation(alggrowth(S,2)));
B := VectorSpace(K,Concatenation(alggrowth(S,6)));
LP := function(i,j)
  return VectorSpace(K,List(Cartesian(Basis(i),Basis(j)),p->p[1]*p[2]-p[2]*p[1]));
end;
LS := function(b,v)
  local L;
  L := [b];
  while Size(L)<3 and not IsTrivial(L[Size(L)]) do
    Add(L,Intersection(b,gengrowth(S,Basis(LP(L[Size(L)],v)),3)));
  od;
  return L;
end;
