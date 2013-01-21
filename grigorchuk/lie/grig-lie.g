
InstallMethod( NormalizedElementOfMagmaRingModuloRelations,
    "for family of free Lie algebra elements, and list",
    true,
    [ IsFamilyElementOfFreeLieAlgebra, IsList ], 0,
    function( Fam, descr )

        local i;

for i in [1,3..Length(descr[2])-1] do
   descr[2][i]:= ExtRepOfObj( descr[2][i] );
od;


return ObjByExtRep( Fam, descr );

    end );

L := FreeLieAlgebra(GF(2),"a","b","c","d");
a := L.1; b := L.2; c := L.3; d := L.4;
quo := function(R,n)
x := a*b;
ox := x*d*a;
oox := ox*b*a*c*a;
ooox := oox*d*a*c*a*b*a*c*a;
oooox := ooox*c*a*c*a*b*a*c*a*d*a*c*a*b*a*c*a;
r := Concatenation([b+c+d,
  b*d,
  x*(a-b),a*d*a,a*d*d,
  x*a*b,
  ox*d*a*b,ox*d*a*d,ox*b*a*d,
  oox*d*a*d,oox*c*a*d,
  oox*d*a*b*a*c,oox*c*a*c*a*c,oox*c*a*c*a*b,
  ooox*b*a*d,
  ooox*b*a*c*a*c,
  ooox*b*a*c*a*b*a*c*a*c,ooox*b*a*c*a*b*a*c*a*b,ooox*c*a*c*a*b*a*c*a*b
],R);
Q := NilpotentQuotientOfFpLieAlgebra(L/r,n);
A := GeneratorsOfAlgebra(Source(Q))[1]^Q;
B := GeneratorsOfAlgebra(Source(Q))[2]^Q;
C := GeneratorsOfAlgebra(Source(Q))[3]^Q;
D := GeneratorsOfAlgebra(Source(Q))[4]^Q;
return Range(Q);
end;
