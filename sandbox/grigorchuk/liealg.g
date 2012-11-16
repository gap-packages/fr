M := JenningsLieAlgebra(GF(2),AsLpGroup(GrigorchukEvilTwin));
a := GeneratorsOfAlgebra(M)[1];
b := GeneratorsOfAlgebra(M)[2];
c := GeneratorsOfAlgebra(M)[3];
d := GeneratorsOfAlgebra(M)[4];

xxx := function() Add(B,[B[Length(B)][1]*a,B[Length(B)][2]*a]); Append(B[Length(B)],List(B[Length(B)/2]{[1..2]},x->x^2)); end;

B := [[a,b,c,d],[a*b,a*c,a*d,c*d]];
Add(B,[B[Length(B)][1]*c,B[Length(B)][2]*d]); xxx();
Add(B,[B[Length(B)][1]*d,B[Length(B)][2]*c]); xxx();
Add(B,[B[Length(B)][1]*c,B[Length(B)][2]*c]); xxx();
Add(B,[B[Length(B)][1]*d,B[Length(B)][2]*d]); xxx();
Add(B,[B[Length(B)][1]*c,B[Length(B)][2]*c]); xxx();
Add(B,[B[Length(B)][1]*d,B[Length(B)][2]*d]); xxx();
Add(B,[B[Length(B)][1]*c,B[Length(B)][2]*c]); xxx();
Add(B,[B[Length(B)][1]*c,B[Length(B)][2]*c]);
B := List([1..Length(B)],i->Basis(Grading(M).hom_components(i),B[i]));
################################################################

L<a,b,c,d> := FreeLieAlgebra(GF(2),4);
bcd := b+c+d;
R := [b*bcd,c*bcd,
      a*b*b,a*c*c,a*d*d,a*b*a,a*c*a,a*d*a,a*bcd*bcd,a*bcd*d,b*c*d,b*c*c,
      a*b*c*a*(b+c),a*c*d*a*d,a*b*c*a*c-a*c*d*a*c,
      a*c*d*a*c*a*d,
      a*b*c*a*d*a*c*a*b,a*b*c*a*d*a*c*a*bcd,a*c*d*a*c*a*c*a*c,
      a*c*d*a*c*a*c*a*d*a*c*a*c,
      a*b*c*a*d*a*c*a*d*a*c*a*d*a*c*a*d,a*b*c*a*d*a*c*a*d*a*c*a*d*a*c*a*bcd,a*c*d*a*c*a*c*a*d*a*c*a*d*a*c*a*(c+d),
      a*c*d*a*c*a*c*a*d*a*c*a*d*a*c*a*d*a*c*a*d*a*c*a*(c+d)];
 time K, B, G, f := NilpotentQuotient(R,33); G;
 