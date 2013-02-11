#############################################################################
##
#W bisets.gi                                                Laurent Bartholdi
##
#Y Copyright (C) 2012-2013, Laurent Bartholdi
##
#############################################################################
##
##  This file implements general bisets.
##
#############################################################################

################################################################
# create new bisets
################################################################
BindGlobal("FRBISET_NEWTYPE@", function()
    return NewType(FRBISET_FAMILY,IsFRBiset and IsFRBisetByFRMachine);
end);

InstallMethod(BisetByFRMachine, [IsFRMachine],
        function(M)
    local b, g;
    b := Objectify(FRBISET_NEWTYPE@(), rec(machine := M));
    if IsInvertible(M) then
        g := SCGroup(M);
    else
        g := SCMonoid(M);
    fi;
    SetLeftActingDomain(b,g);
    SetRightActingDomain(b,g);
    return b;
end);

InstallMethod(BisetByFRSemigroup, [IsFRGroup],
        function(G)
    local b;
    b := Objectify(FRBISET_NEWTYPE@(), rec(machine := UnderlyingFRMachine(G)));
    SetLeftActingDomain(b,G);
    SetRightActingDomain(b,G);
    return b;
end);

InstallMethod(ViewString, "(FR) for a biset",
        [IsFRBiset],
        function(b)
    local s, v;
    s := "<";
    for v in [[IsLeftFree,"left-free"],
            [IsRightFree,"right-free"],
            [IsLeftTransitive,"left-transitive"],
            [IsRightTransitive,"right-transitive"]] do
        if Tester(v[1])(b) then
            if not v[1](b) then Append(s,"non-"); fi;
            Append(s,v[2]);
        fi;
    od;
    if not s[Length(s)] in " <" then Append(s," "); fi;
    if HasName(LeftActingDomain(b)) and HasName(RightActingDomain(b)) then
        APPEND@(s,Name(LeftActingDomain(b)),"-",Name(RightActingDomain(b)),"-");
    fi;
    Append(s,"biset>");
    return s;
end);
InstallMethod(DisplayString, "(FR) for a biset",
        [IsFRBiset and IsFRBisetByFRMachine],
        function(b)
    return Concatenation("Biset:\n",DisplayString(b!.machine));
end);
INSTALLPRINTERS@(IsFRBiset);

################################################################
# biset elements
################################################################
InstallMethod(BisetElement, "(FR) for a machine biset",
        [IsFRBisetByFRMachine,IsMultiplicativeElement,IsObject],
        function(biset,g,x)
    local type, elt;
    
    type := NewType(FRBISET_FAMILY,IsBisetElement and IsBisetElementByPair);
    elt := rec(biset := biset, groupelement := g, letter := x);
    return Objectify(type,elt);
end);

InstallMethod(BisetElement, "(FR) for a homomorphism biset",
        [IsFRBisetByHomomorphism,IsMultiplicativeElement],
        function(biset,g)
    local type, elt;
    
    type := NewType(FRBISET_FAMILY,IsBisetElement and IsBisetElementByElement);
    elt := rec(biset := biset, groupelement := g);
    return Objectify(type,elt);
end);

InstallMethod(EQ, "(FR) for two biset elements",
        [IsBisetElementByPair,IsBisetElementByPair],
        function(e,f)
    return e!.biset=f!.biset and e!.groupelement=f!.groupelement and e!.letter=f!.letter;
end);
InstallMethod(LT, "(FR) for two biset elements",
        [IsBisetElementByPair,IsBisetElementByPair],
        function(e,f)
    return e!.biset<f!.biset or (e!.biset=f!.biset and (e!.letter<f!.letter or (e!.letter=f!.letter and e!.groupelement<f!.groupelement)));
end);
         
InstallMethod(EQ, "(FR) for two biset elements",
        [IsBisetElementByElement,IsBisetElementByElement],
        function(e,f)
    return e!.biset=f!.biset and e!.groupelement=f!.groupelement;
end);
InstallMethod(LT, "(FR) for two biset elements",
        [IsBisetElementByElement,IsBisetElementByElement],
        function(e,f)
    return e!.biset<f!.biset or (e!.biset=f!.biset and e!.groupelement<f!.groupelement);
end);
         
InstallMethod(ViewString, "(FR) for a biset element by pair",
        [IsBisetElement and IsBisetElementByPair],
        function(e)
    return CONCAT@FR("(",e!.groupelement,"*",e!.letter,")");
end);
InstallMethod(PrintString, "(FR) for a biset element by pair",
        [IsBisetElement and IsBisetElementByPair],
        function(e)
    return CONCAT@FR("BisetElement(",e!.biset,",",e!.groupelement,",",e!.letter,")");
end);
INSTALLPRINTERS@(IsBisetElement);

################################################################
# bases
################################################################
InstallMethod(Basis, "(FR) for a machine biset",
        [IsFRBiset],
        LeftBasis);

InstallMethod(LeftBasis, "(FR) for a machine biset",
        [IsFRBisetByFRMachine],
        CanonicalBasis);

InstallMethod(CanonicalBasis, "(FR) for a machine biset",
        [IsFRBisetByFRMachine],
        function(b)
    return Objectify(NewType(FRBISET_FAMILY,IsLeftBisetBasis and IsCanonicalBasis),
                     List(AlphabetOfFRObject(b!.machine),i->BisetElement(b,One(LeftActingDomain(b)),i)));
end);

InstallMethod(ELM_LIST, "(FR) for a biset basis",
        [IsBisetBasis,IsPosInt],
        function(basis,i)
    return basis![i];
end);

InstallMethod(BasisVectors, "(FR) for a biset basis",
        [IsBisetBasis],
        function(basis)
    local v, i;
    v := [];
    i := 1;
    while IsBound(basis![i]) do
        Add(v,basis![i]);
        i := i+1;
    od;
    return v;
end);

InstallMethod(ViewString, "(FR) for a biset basis",
        [IsLeftBisetBasis],
        function(basis)
    return "LeftBasis(...)";
end);

InstallMethod(ViewString, "(FR) for a biset basis",
        [IsRightBisetBasis],
        function(basis)
    return "RightBasis(...)";
end);
INSTALLPRINTERS@(IsBisetBasis);

#E bisets.gi. . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
