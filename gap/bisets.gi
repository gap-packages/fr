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
InstallMethod(BisetByFRMachine, [IsFRMachine],
        function(M)
    local b, g;
    
    if IsMealyMachine(M) then
        if IsInvertible(M) then
            M := AsGroupFRMachine(M);
        else
            M := AsMonoidFRMachine(M);
        fi;
    fi;
    b := Objectify(NewType(FRBISET_FAMILY,IsFRBiset and IsFRBisetByFRMachineRep),
                 rec(machine := M));
    g := StateSet(M);
    SetIsLeftFree(b,true);
    SetLeftActingDomain(b,g);
    SetRightActingDomain(b,g);
    return b;
end);

InstallMethod(BisetByFRSemigroup, [IsFRGroup],
        function(G)
    local b;
    b := Objectify(NewType(FRBISET_FAMILY,IsFRBiset and IsFRBisetByFRSemigroupRep),
                 rec(semigroup := G));
    SetIsLeftFree(b,true);
    SetLeftActingDomain(b,G);
    SetRightActingDomain(b,G);
    return b;
end);

InstallMethod(DualBiset, "(FR) for a biset",
        [IsFRBiset],
        function(b)
    Error("not yet done");
end);

InstallMethod(TensorProductOp, "(FR) for a sequence of bisets and a biset",
        [IsList,IsFRBiset],
        function(b,b0)
    Error("not yet done");
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
        [IsFRBiset and IsFRBisetByFRMachineRep],
        function(b)
    return Concatenation("Biset:\n",DisplayString(b!.machine));
end);
INSTALLPRINTERS@(IsFRBiset);

################################################################
# biset elements
################################################################
InstallMethod(BisetElement, "(FR) for a machine biset",
        [IsFRBiset,IsMultiplicativeElement,IsObject],
        function(biset,g,x)
    local type, elt;
    
    type := NewType(FRBISET_FAMILY,IsBisetElement and IsBisetElementByPair);
    elt := rec(biset := biset, element := g, letter := x);
    return Objectify(type,elt);
end);

InstallMethod(BisetElement, "(FR) for a homomorphism biset",
        [IsFRBisetByHomomorphismRep,IsMultiplicativeElement],
        function(biset,g)
    local type, elt;
    
    type := NewType(FRBISET_FAMILY,IsBisetElement and IsBisetElementByElement);
    elt := rec(biset := biset, element := g);
    return Objectify(type,elt);
end);

InstallMethod(\=, "(FR) for two biset elements",
        [IsBisetElementByPair,IsBisetElementByPair],
        function(e,f)
    return e!.biset=f!.biset and e!.element=f!.element and e!.letter=f!.letter;
end);
InstallMethod(\<, "(FR) for two biset elements",
        [IsBisetElementByPair,IsBisetElementByPair],
        function(e,f)
    return e!.biset<f!.biset or (e!.biset=f!.biset and (e!.letter<f!.letter or (e!.letter=f!.letter and e!.element<f!.element)));
end);
         
InstallMethod(\=, "(FR) for two biset elements",
        [IsBisetElementByElement,IsBisetElementByElement],
        function(e,f)
    return e!.biset=f!.biset and e!.element=f!.element;
end);
InstallMethod(\<, "(FR) for two biset elements",
        [IsBisetElementByElement,IsBisetElementByElement],
        function(e,f)
    return e!.biset<f!.biset or (e!.biset=f!.biset and e!.element<f!.element);
end);

InstallOtherMethod(\*, "(FR) for a biset element and a semigroup element",
        [IsBisetElementByPair,IsMultiplicativeElement],
        function(b,g)
    local biset, newelement;
    biset := b!.biset;
    if not g in LeftActingDomain(biset) then TryNextMethod(); fi;
    if IsFRBisetByFRMachineRep(biset) then
        return BisetElement(biset,b!.element*Transition(biset!.machine,g,b!.letter),Output(biset!.machine,g,b!.letter));
    elif IsFRElement(g) then
        return BisetElement(b!.biset,b!.element*State(g,b!.letter),b!.letter^g);
    else
        TryNextMethod();
    fi;
end);
InstallOtherMethod(\*, "(FR) for a semigroup element and a biset element",
        [IsMultiplicativeElement,IsBisetElementByPair],
        function(g,b)
    if not g in RightActingDomain(b!.biset) then TryNextMethod(); fi;
    return BisetElement(b!.biset,g*b!.element,b!.letter);
end);
  
InstallMethod(ViewString, "(FR) for a biset element by pair",
        [IsBisetElement and IsBisetElementByPair],
        function(e)
    return CONCAT@("(",e!.element,"*",e!.letter,")");
end);
InstallMethod(PrintString, "(FR) for a biset element by pair",
        [IsBisetElement and IsBisetElementByPair],
        function(e)
    return CONCAT@("BisetElement(",e!.biset,",",e!.element,",",e!.letter,")");
end);
INSTALLPRINTERS@(IsBisetElement);

################################################################
# bases
################################################################
InstallMethod(Basis, "(FR) for a machine biset",
        [IsFRBiset],
        LeftBasis);

InstallMethod(LeftBasis, "(FR) for a machine biset",
        [IsFRBiset and IsLeftFree],
        CanonicalBasis);

BindGlobal("CANONICALBASIS@", function(b,alphabet)
    return Objectify(NewType(FRBISET_FAMILY,IsLeftBisetBasis and IsCanonicalBasis),
                     List(alphabet,i->BisetElement(b,One(LeftActingDomain(b)),i)));
end);

InstallMethod(CanonicalBasis, "(FR) for a machine biset",
        [IsFRBisetByFRMachineRep],
        b->CANONICALBASIS@(b,AlphabetOfFRObject(b!.machine)));

InstallMethod(CanonicalBasis, "(FR) for a machine biset",
        [IsFRBisetByFRSemigroupRep],
        b->CANONICALBASIS@(b,AlphabetOfFRSemigroup(b!.semigroup)));

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

InstallMethod(WreathRecursion, "(FR) for a biset",
        [IsFRBiset and IsFRBisetByFRMachineRep],
        b->WreathRecursion(b!.machine));
InstallMethod(WreathRecursion, "(FR) for a biset",
        [IsFRBiset],
        b->WreathRecursion(b,Basis(b)));
InstallMethod(WreathRecursion, "(FR) for a biset with basis",
        [IsFRBiset,IsLeftBisetBasis],
        function(biset,basis)
    Error("not yet done");
end);

InstallMethod(FRMachineOfBiset, "(FR) for a machine biset",
        [IsFRBiset and IsFRBisetByFRMachineRep],
        b->b!.machine);
InstallMethod(FRMachineOfBiset, "(FR) for a biset",
        [IsFRBiset],
        b->FRMachine(b,Basis(b)));
InstallMethod(FRMachine, "(FR) for a biset",
        [IsFRBiset],
        FRMachineOfBiset);
InstallMethod(FRMachine, "(FR) for a biset with basis",
        [IsFRBiset,IsLeftBisetBasis],
        function(biset,basis)
    Error("not yet done");
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

