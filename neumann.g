G := Group((1,2,3),(4,5,6));

pi := EpimorphismSchurCover(G);
GT := Source(pi); # Schur cover
H2 := Kernel(pi); # element in H2 is homology class.
H2s := IrreducibleRepresentations(H2); # element in H2s is cohomology class.

DeclareProperty("IsRepresentation",IsGroupHomomorphism);
DeclareAttribute("Contragredient",IsRepresentation);
Perform(H2s,function(r) SetIsRepresentation(r,true); end);
InstallOtherMethod(TensorProductOp,[IsList,IsRepresentation],
        function(list,tag)
    local r, s, t, g, m, phi;
    s := Source(tag);
    if not ForAll(list,r->Source(r)=s) then
        TryNextMethod();
    fi;
    r := [];
    for g in GeneratorsOfGroup(s) do
        m := [[1]];
        for phi in list do
            m := KroneckerProduct(m,Image(phi,g));
        od;
        Add(r,m);
    od;
    t := GroupHomomorphismByImages(s,Group(r),GeneratorsOfGroup(s),r);
    SetIsRepresentation(t,true);
    return t;
end);
InstallMethod(Contragredient,[IsRepresentation],
        function(phi)
    local gens;
    gens := GeneratorsOfGroup(Source(phi));
    return GroupHomomorphismByImages(Source(phi),Range(phi),gens,List(gens,g->Inverse(TransposedMat(Image(phi,g)))));
end);
    
canonical := function(g,h)
    # tautological cocycle GxG -> H2
    return PreImagesRepresentative(pi,g)*PreImagesRepresentative(pi,h)/PreImagesRepresentative(pi,g*h);
end;
Pcgs(H2)!.rep := [];
for g in Pcgs(H2) do
    Add(Pcgs(H2)!.rep,First(Cartesian(G,G),p->canonical(p[1],p[2])=g));
od;

repcocycle := function(chi)
    # takes an element in H2*, returns a function GxG -> C*
    return function(g,h)
        return Image(chi,canonical(g,h));
    end;
end;

cocyclerep := function(f)
    # takes a function GxG -> C*, constructs an element of H2*
    local v;
    v := List(Pcgs(H2)!.rep,p->f(p[1],p[2]));
    return GroupHomomorphismByImages(H2,Group(v),Pcgs(H2),v);
end;
