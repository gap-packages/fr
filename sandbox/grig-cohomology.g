LoadPackage("happrime");
N := 8;
q := [];
g := PcGroup(GrigorchukGroup,6);
q[N] := g/LowerCentralSeries(g)[N];
f := [];
for i in [N-1,N-2..1] do
    f[i+1] := NaturalHomomorphismByNormalSubgroup(q[i+1],LowerCentralSeries(q[i+1])[i]);
    q[i] := Image(f[i+1]);
od;
r := List([1..N],i->ResolutionNilpotentGroup(q[i],4));
Z_r := List(r,TensorWithIntegers);
Z2_r := List(r,x->TensorWithIntegersModP(x,2));
ZG_f := List([1..N-1],i->EquivariantChainMap(r[i+1],r[i],f[i+1]));
Z_f := List(ZG_f,TensorWithIntegers);
Z2_f := List(ZG_f,f->TensorWithIntegersModP(f,2));

#List([1..5],i->Collected(Homology(Z_r[3],i)));
#List([1..5],i->Collected(AbelianInvariants(Image(Homology(Z_f[2],i)))));
#List([2..N-1],i->Collected(AbelianInvariants(Image(Homology(Z2_f[i],2))))[1][2]);
