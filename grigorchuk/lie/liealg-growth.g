n := 12;
L := FreeLieAlgebra(GF(2),4);
a := L.1;
b := L.2;
c := L.3;
d := L.4;
L := L / [b+c+d,b*c,a*d*d,a*b*b,a*b*a,a*d*a,
          a*b*d*a*d,a*b*d*a*b*a*d,a*b*d*a*b*a*b*a*c];
Q := NilpotentQuotientOfFpLieAlgebra(L,n);
a := GeneratorsOfAlgebra(L)[1]^Q;
b := GeneratorsOfAlgebra(L)[2]^Q;
c := GeneratorsOfAlgebra(L)[3]^Q;
d := GeneratorsOfAlgebra(L)[4]^Q;
S := LieLowerCentralSeries(Range(Q));
Print(List([1..n],i->Dimension(S[i])-Dimension(S[i+1])),"\n");

M := JenningsLieAlgebra(GF(2),AsLpGroup(GrigorchukGroup));
a := GeneratorsOfAlgebra(M)[1];
b := GeneratorsOfAlgebra(M)[2];
c := GeneratorsOfAlgebra(M)[3];
d := GeneratorsOfAlgebra(M)[4];

################################################################

L<a,b,c> := FreeLieAlgebra(GF(2),3);
d := b+c;
R := [b+c+d,b*c,a*b*a,a*d*a,b*a*b,
      d                               *                               a                               *                               d,
      d               *               a               *               c               *               a               *               d,
      d       *       a       *       c       *       a       *       b       *       a       *       c       *       a       *       c,
      d   *   a   *   c   *   a   *   b   *   a   *   c   *   a   *   d   *   a   *   c   *   a   *   b   *   a   *   c   *   a   *   b,
      d * a * c * a * b * a * c * a * d * a * c * a * b * a * c * a * c * a * c * a * b * a * c * a * d * a * c * a * b * a * c * a * d,
      d*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a*c*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a*b*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a*c*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a*c,
      d       *       a       *       c       *       a       *       c       *       a       *       d,
      d   *   a   *   c   *   a   *   b   *   a   *   c   *   a   *   b   *   a   *   c   *   a   *   c,
      d * a * c * a * b * a * c * a * d * a * c * a * b * a * c * a * d * a * c * a * b * a * c * a * b,
      d*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a*c*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a*b*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a*d];
time K, B, G, f := NilpotentQuotient(R,98); G;
