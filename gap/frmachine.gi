#############################################################################
##
#W frmachine.gi                                             Laurent Bartholdi
##
#Y Copyright (C) 2006-2013, Laurent Bartholdi
##
#############################################################################
##
##  This file implements the category of functionally recursive machines.
##
#############################################################################

############################################################################
##
#O AlphabetOfFRObject
#V FR_FAMILIES
#O FRMFamily(alphabet)
##
InstallMethod(AlphabetOfFRObject, "(FR) for an FR object",
        [IsFRObject],
        function(M)
    local a;
    a := FamilyObj(M)!.alphabet;
    IsRange(a);
    return a;
end);

INSTALLPRINTERS@(IsFRObject);

InstallMethod(FRMFamily, "(FR) for an alphabet",
        [IsListOrCollection],
        function(d)
    local i;
    for i in FR_FAMILIES do
        if i[1] = d then return i[2]; fi;
    od;
    i := NewFamily(Concatenation("FRMachine(",String(d),")"), IsFRMachine);
    i!.standard := Size(d)<2^28 and d=[1..Size(d)];
    if i!.standard then
        i!.alphabet := [1..Size(d)];
    else
        i!.alphabet := d;
        i!.a2n := x->Position(Enumerator(d),x);
        i!.n2a := x->Enumerator(d)[x];
    fi;
    ConvertToRangeRep(i!.alphabet);
    MakeImmutable(i!.alphabet);
    Add(FR_FAMILIES,[d,i]);
    return i;
end);

InstallMethod(IsGroupFRMachine, [IsFRMachine], ReturnFalse);
InstallMethod(IsMonoidFRMachine, [IsFRMachine], ReturnFalse);
InstallMethod(IsSemigroupFRMachine, [IsFRMachine], ReturnFalse);
#############################################################################

#############################################################################
##
#O FRMachine(Transitions, Output)
#O FRMachine(Names, Transitions, Output)
#O FRMachineNC(Family, [Semi]Group, Transitions, Output)
#O FRMachine([Semi]Group, Transitions, Output)
##
InstallOtherMethod(FRMachineNC, "(FR) for a family, a free group, a list of transitions and a list of outputs",
        [IsFamily, IsGroup, IsList, IsList],
        function(fam,free,transitions,output)
    local M, F;
    F := FamilyObj(Representative(free));
    M := Objectify(NewType(fam, IsGroupFRMachine and IsFRMachineStdRep),
                 rec(free := free,
                     pack := l->AssocWordByLetterRep(F,l),
                     transitions := Immutable(transitions),
                     output := Immutable(output)));
    SetIsInvertible(M, true);
    return M;
end);

InstallOtherMethod(FRMachineNC, "(FR) for a family, a free semigroup, a list of transitions and a list of outputs",
        [IsFamily, IsFreeSemigroup, IsList, IsList],
        function(fam,free,transitions,output)
    local M, F;
    F := FamilyObj(Representative(free));
    M := Objectify(NewType(fam, IsSemigroupFRMachine and IsFRMachineStdRep),
                 rec(free := free,
                     pack := l->AssocWordByLetterRep(F,l),
                     transitions := Immutable(transitions),
                     output := Immutable(output)));
    return M;
end);

InstallOtherMethod(FRMachineNC, "(FR) for a family, a free monoid, a list of transitions and a list of outputs",
        [IsFamily, IsFreeMonoid, IsList, IsList],
        function(fam,free,transitions,output)
    local M, F, T;
    F := FamilyObj(Representative(free));
    M := Objectify(NewType(fam, IsMonoidFRMachine and IsFRMachineStdRep),
                 rec(free := free,
                     pack := l->AssocWordByLetterRep(F,l),
                     transitions := Immutable(transitions),
                     output := Immutable(output)));
    return M;
end);

BindGlobal("COPYFRMACHINE@", function(m)
    return Objectify(NewType(FamilyObj(m), First([IsGroupFRMachine,IsMonoidFRMachine,IsSemigroupFRMachine],p->Tester(p)(m) and p(m)) and IsFRMachineStdRep),
                   rec(free := m!.free,
                       pack := m!.pack,
                       transitions := m!.transitions,
                       output := m!.output));
end);

BindGlobal("ANY2OUT@", function(x,n)
    if IsList(x) then
        return x;
    elif IsTransformation(x) then
        return ListTransformation(x,n);
    elif IsPerm(x) then
        return ListPerm(x,n);
    fi;
end);

BindGlobal("CHECKLENGTHSCONTENTS@", function(t, transitions, output)
    # check validity of arguments;
    # unpack FR elements contained in the transitions;
    # set t.F
    local i, j, k, x, e;
    if Length(transitions)<>Length(output) then
        Error("<Transitions> and <Output> must have the same length\n");
        return fail;
    fi;
    if not ForAll(transitions, IsList) or
       ForAny(transitions, r->Length(r)<>Length(transitions[1])) then
        Error("All rows of <Transitions> must be lists of the same length\n");
        return fail;
    fi;
    t.F := FRMFamily([1..Length(transitions[1])]);
    t.transitions := StructuralCopy(transitions);
    t.output := ShallowCopy(output);
    for x in t.transitions do for x in x do if IsList(x) then
        i := 1; while i <= Length(x) do
            if IsFRElement(x[i]) then
                if IsMealyElement(x[i]) then
                    e := AsSemigroupFRElement(x[i]);
                else
                    e := x[i];
                fi;
                k := Length(t.transitions);
                for j in UnderlyingFRMachine(e)!.transitions do
                    Add(t.transitions,List(j,w->List(LetterRepAssocWord(w),i->i+SignInt(i)*k)));
                od;
                Append(t.output, UnderlyingFRMachine(e)!.output);
                Remove(x,i); i := i-1;
                for j in LetterRepAssocWord(InitialState(e)) do
                    i := i+1;
                    Add(x,j+SignInt(j)*k,i);
                od;
            elif IsInt(x[i]) and AbsInt(x[i]) in [1..Length(t.transitions)] then;
            else
                Error("Entry ",i," of <Transitions> is not in the state set\n");
                return fail;
            fi;
            i := i+1;
        od;
    elif not IsAssocWord(x) then
        Error("Transitions must be associative words or lists");
    fi; od; od;

    # clean up t.output, set t.invertible
    t.invertible := true;
    for i in [1..Length(t.output)] do
        t.output[i] := ANY2OUT@(t.output[i],Length(t.F!.alphabet));
        if Set(t.output[i])<>t.F!.alphabet then
            t.invertible := false;
        fi;
    od;
    for i in t.output do
        if not IsSubset(t.F!.alphabet,i) then
            Error("Entry ",i," of <Output> is not in alphabet ",t.F!.alphabet,"\n");
            return fail;
        fi;
    od;
end);

InstallMethod(FRMachine, "(FR) for a list of transitions and a list of outputs",
        [IsList, IsList],
        function(transitions, output)
    local G, elG, t;
    t := rec();
    CHECKLENGTHSCONTENTS@(t, transitions, output);
    if t.invertible then
        G := FreeGroup(Length(t.transitions));
    else
        G := FreeMonoid(Length(t.transitions));
    fi;
    elG := FamilyObj(Representative(G));
    return FRMachineNC(t.F, G, List(t.transitions,t->List(t,w->AssocWordByLetterRep(elG,w))),t.output);
end);

InstallMethod(FRMachine, "(FR) for a list of names, a list of transitions and a list of outputs",
        [IsList, IsList, IsList],
        function(names, transitions, output)
    local G, elG, t, n;
    t := rec();
    if not ForAll(names,IsString) then
        Error("<names> should be a list of strings, and not ", names,"\n");
    fi;
    CHECKLENGTHSCONTENTS@(t, transitions, output);
    if Length(names)>Length(t.transitions) then
        Error("Too many names supplied to FRMachine()\n");
    elif Length(names)<Length(t.transitions) then
        n := Concatenation(names,List([1..Length(t.transitions)-Length(names)],i->Concatenation("__",String(i))));
    else
        n := names;
    fi;
    if t.invertible then
        G := FreeGroup(n);
    else
        G := FreeMonoid(n);
    fi;
    elG := FamilyObj(Representative(G));
    return FRMachineNC(t.F, G, List(t.transitions,t->List(t,w->AssocWordByLetterRep(elG,w))),t.output);
end);

InstallMethod(FRMachine, "(FR) for a free [semi]group, a list of transitions and a list of outputs",
        [IsSemigroup, IsList, IsList],
        function(free,transitions,output)
    local t, elfree, r, i;
    t := rec();
    CHECKLENGTHSCONTENTS@(t, transitions, output);
    if IsGroup(free) and not t.invertible then
        Error("Outputs must be invertible in group FR machine: ",t.output);
    fi;
    elfree := FamilyObj(Representative(free));
    for r in t.transitions do for i in [1..Length(r)] do
        if IsList(r[i]) then r[i] := AssocWordByLetterRep(elfree,r[i]); fi;
    od; od;
    r := FRMachineNC(t.F,free,t.transitions,t.output);
    if Length(t.transitions)<>Length(GeneratorsOfFRMachine(r)) then
        Error("<Transition> and <Output> should have same length as ",free,"'s rank\n");
    fi;
    return r;
end);

#############################################################################
##
#A  GeneratorsOfFRMachine(FRMachine)
##
InstallMethod(GeneratorsOfFRMachine, "(FR) for a group FR machine",
        [IsGroupFRMachine],
        M->GeneratorsOfGroup(M!.free));

InstallMethod(GeneratorsOfFRMachine, "(FR) for a semigroup FR machine",
        [IsSemigroupFRMachine],
        M->GeneratorsOfSemigroup(M!.free));

InstallMethod(GeneratorsOfFRMachine, "(FR) for a monoid FR machine",
        [IsMonoidFRMachine],
        M->GeneratorsOfMonoid(M!.free));

InstallMethod(StateSet, "(FR) for an FR machine",
        [IsFRMachine and IsFRMachineStdRep],
        M->M!.free);
#############################################################################

#############################################################################
##
#M  ViewObj(FRMachine)
#M  String(FRMachine)
#M  Display(FRMachine)
##
InstallMethod(ViewString, "(FR) for an FR machine",
        [IsFRMachine and IsFRMachineStdRep],
        function(M)
    return CONCAT@("<FR machine with alphabet ", AlphabetOfFRObject(M), " on ", StateSet(M), ">");
end);

InstallMethod(String, "(FR) for an FR machine",
        [IsFRMachine and IsFRMachineStdRep],
        function(M)
    return CONCAT@("FRMachine(...,",M!.output,")");
end);

BindGlobal("DISPLAYFRMACHINE@", function(M)
    local a, i, j, g, alen, slen, glen, ablank, sblank, gblank, arule, grule, srule, StringId, s;
    a := AlphabetOfFRObject(M);
    g := GeneratorsOfFRMachine(M);
    s := "";
    StringId := function(arg)
        local s;
        s := CallFuncList(String,arg);
        if s="<identity ...>" then
            s := "<id>"; if Length(arg)=2 then s := String(s,arg[2]); fi;
        fi;
        return s;
    end;
    alen := LogInt(Maximum(a),10)+3;
    ablank := ListWithIdenticalEntries(alen,' ');
    arule := ListWithIdenticalEntries(alen,'-');
    if g=[] then
        glen := 2;
        slen := List(a,i->1);
    else
        glen := Maximum(List(g,t->Length(StringId(t))))+1;
        slen := List(a,i->Maximum(List(g,t->Length(StringId(Transition(M,t,i)))))+1);
    fi;
    gblank := ListWithIdenticalEntries(glen,' ');
    grule := ListWithIdenticalEntries(glen,'-');
    sblank := List(a,i->ListWithIdenticalEntries(slen[i],' '));
    srule := List(a,i->ListWithIdenticalEntries(slen[i],'-'));

    if IsGroupFRMachine(M) then
        s := " G";
    elif IsMonoidFRMachine(M) then
        s := " M";
    else s := " S"; fi;
    APPEND@(s,gblank{[3..glen]}," |");
    for i in [1..Length(a)] do APPEND@(s,sblank[i],String(a[i],-alen)," "); od;
    APPEND@(s,"\n");
    APPEND@(s,grule,"-+");
    for i in [1..Length(a)] do APPEND@(s,srule[i],arule,"+"); od;
    APPEND@(s,"\n");
    for i in [1..Length(g)] do
        APPEND@(s,StringId(g[i],glen)," |");
        for j in [1..Length(a)] do
            APPEND@(s,StringId(M!.transitions[i][j],slen[j]),",",String(M!.output[i][j],-alen));
        od;
        APPEND@(s,"\n");
    od;
    APPEND@(s,grule,"-+");
    for i in [1..Length(a)] do APPEND@(s,srule[i],arule,"+"); od;
    APPEND@(s,"\n");
    return s;
end);

InstallMethod(DisplayString, "(FR) for an FR machine",
        [IsFRMachine and IsFRMachineStdRep],
        DISPLAYFRMACHINE@);

INSTALLPRINTERS@(IsFRMachine);
#############################################################################
##
#M One(FRMachine)
##
InstallMethod(OneOp, "(FR) for an FR machine",
        [IsFRMachine],
        M->FRMachine([ListWithIdenticalEntries(Size(AlphabetOfFRObject(M)),[1])],[()]));
#############################################################################

#############################################################################
##
#M Zero(FRMachine)
##
InstallOtherMethod(ZeroOp, "(FR) for an FR machine",
        [IsFRMachine],
        M->FRMachineNC(FamilyObj(M),FreeGroup(0),[],[]));
#############################################################################

#############################################################################
##
#M InverseOp(FRMachine)
##
InstallTrueMethod(IsInvertible, IsGroupFRMachine);

BindGlobal("ISINVERTIBLE@", function(l)
    return Set(l)=[1..Length(l)];
end);

BindGlobal("INVERSE@", function(l) # inverse of transformation, given as list
    local r;
    r := [];
    r{l} := [1..Length(l)];
    return r;
end);

BindGlobal("ISONE@", function(l) # identity mapping, given as list
    return l=[1..Length(l)];
end);

BindGlobal("PREIMAGE@", Position); # preimage of point under transformation

InstallMethod(InverseOp, "(FR) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    local N;
    N := FRMachineNC(FamilyObj(M), M!.free,
                 List([1..Length(M!.transitions)], i->M!.transitions[i]{M!.output[i]}),
                 List(M!.output, INVERSE@));
    SetInverse(M,N);
    SetInverse(N,M);
    return N;
end);

InstallMethod(IsReversible, "(FR) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    local a, hom;
    for a in AlphabetOfFRObject(M) do
        hom := GroupHomomorphismByImages(StateSet(M),StateSet(M),GeneratorsOfFRMachine(M),M!.transitions{[1..Length(M!.transitions)]}[a]);
        if not IsBijective(hom) then
            return false;
        fi;
    od;
    return true;
end);
#############################################################################

#############################################################################
##
#M Products
##
BindGlobal("SET_NAME@", function(args,sep,obj)
    local i, s, n;
    for i in args do if not HasName(i) then return; fi; od;
    s := ShallowCopy(Name(args[1]));
    for i in [2..Length(args)] do Append(s,sep); Append(s,Name(args[i])); od;
    SetName(obj,s);
end);

BindGlobal("MAKENAMESUNIQUE@", function(sgen)
    local i, j, nonunique;
    nonunique := Set(Filtered(Collected(Concatenation(sgen)),x->x[2]>1),x->x[1]);
    RemoveSet(nonunique,"<identity ...>");
    for i in [1..Length(sgen)] do
        for j in [1..Length(sgen[i])] do
            if sgen[i][j] in nonunique then
                sgen[i][j] := Concatenation(sgen[i][j],".",String(i));
            fi;
        od;
    od;
end);

BindGlobal("LARGESTDENOMINATOR@", function(arg)
    # returns homomorphisms from all its arguments' free stateset to
    # a free object of highest structure (group > monoid > semigroup).
    # the last entry in the returned list is a list of appropriate generators
    # for each argument's free object.

    local c, d, f, i, iso, states, gen, sgen, subgen, shift;

    d := Length(arg);
    states := List(arg,x->x!.free);
    c := List(states,IdentityMapping);
    gen := List(arg,GeneratorsOfFRMachine);
    sgen := List(gen,x->List(x,String));

    if ForAll(states,IsGroup) then
        MAKENAMESUNIQUE@(sgen);
        f := FreeGroup(Concatenation(sgen));
        c := [];
        shift := 0;
        for i in [1..d] do
            Add(c,GroupHomomorphismByImages(states[i],f,
                    gen[i],GeneratorsOfGroup(f){shift+[1..Length(gen[i])]}));
            shift := shift + Length(gen[i]);
        od;
    elif ForAll(states,IsMonoid) then
        c := List(states,IdentityMapping);
        for i in [1..d] do
            if IsGroup(states[i]) then
                c[i] := IsomorphismFpMonoidInversesFirst(states[i]);
                gen[i] := GeneratorsOfMonoid(states[i]);
                sgen[i] := List(GeneratorsOfMonoid(Range(IsomorphismFpMonoidInversesFirst(FreeGroup(sgen[i])))),String);
                c[i] := c[i]*MappingByFunction(Range(c[i]),FreeMonoidOfFpMonoid(Range(c[i])),UnderlyingElement);
            fi;
        od;
        MAKENAMESUNIQUE@(sgen);
        f := FreeMonoid(Concatenation(sgen));
        shift := 0;
        for i in [1..d] do
            c[i] := c[i]*MagmaHomomorphismByImagesNC(Range(c[i]),f,
                            GeneratorsOfMonoid(f){shift+[1..Length(gen[i])]});
            shift := shift + Length(gen[i]);
        od;
    else
        c := List(states,IdentityMapping);
        for i in [1..d] do
            if IsGroup(states[i]) then
                c[i] := IsomorphismFpSemigroup(states[i]);
                gen[i] := GeneratorsOfSemigroup(states[i]);
                sgen[i] := List(GeneratorsOfSemigroup(Range(IsomorphismFpSemigroup(FreeGroup(sgen[i])))),String);
                c[i] := c[i]*MappingByFunction(Range(c[i]),FreeSemigroupOfFpSemigroup(Range(c[i])),UnderlyingElement);
            elif IsMonoid(states[i]) then
                gen[i] := GeneratorsOfSemigroup(states[i]);
                sgen[i] := Concatenation(["<identity ...>"],sgen[i]);
                c[i] := states[i]/[];
                iso := IsomorphismFpSemigroup(c[i]);
                c[i] := NaturalHomomorphismByGenerators(states[i],c[i])*iso*MappingByFunction(Range(iso),FreeSemigroupOfFpSemigroup(Range(iso)),UnderlyingElement);
            fi;
        od;
        MAKENAMESUNIQUE@(sgen);
        f := FreeSemigroup(Concatenation(sgen));
        shift := 0;
        for i in [1..d] do
            c[i] := c[i]*MagmaHomomorphismByImagesNC(Range(c[i]),f,
                            GeneratorsOfSemigroup(f){shift+[1..Length(gen[i])]});
            shift := shift + Length(gen[i]);
        od;
    fi;
    Add(c,gen);
    return c;
end);

BindGlobal("FRMSUM@", function(arg)
    local c, gen, trans, out, i, j, sum;

    c := CallFuncList(LARGESTDENOMINATOR@,arg);
    gen := Remove(c);

    trans := [];
    out := [];
    for i in [1..Length(arg)] do
        for j in gen[i] do
            Add(trans,List(AlphabetOfFRObject(arg[i]),a->Transition(arg[i],j,a)^c[i]));
            Add(out,Output(arg[i],j));
        od;
    od;
    sum := FRMachineNC(FamilyObj(arg[1]),Range(c[1]),trans,out);
    SetCorrespondence(sum,c);
    SET_NAME@(arg,"+",sum);
    return sum;
end);

BindGlobal("FRMMINSUM@", function(left,right)
    local sum, r;
    sum := FRMSUM@(left,right);
    r := Minimized(sum);
    r!.Correspondence := List(Correspondence(sum),x->x*Correspondence(r));
    return r;
end);

InstallMethod(\+, "(FR) for two FR machines",
        IsIdenticalObj,
        [IsFRMachine and IsFRMachineStdRep, IsFRMachine and IsFRMachineStdRep],
        FRMSUM@);

InstallMethod(\*, "(FR) for two FR machines",
        IsIdenticalObj,
        [IsFRMachine and IsFRMachineStdRep, IsFRMachine and IsFRMachineStdRep],
        FRMSUM@);

InstallMethod(TensorSumOp, "(FR) for two FR machines",
        [IsList, IsFRMachine and IsFRMachineStdRep],
        function(M, N)
    local trans, out, t, o, i, j, x, d, s;

    while ForAny(M,x->x!.free<>N!.free) do
        Error("All machines should have same stateset");
    od;

    trans := [];
    out := [];
    trans := [];
    out := [];
    for i in [1..Length(N!.output)] do
        t := [];
        o := [];
        d := 0;
        for j in [1..Length(M)] do
            Append(t,M[j]!.transitions[i]);
            Append(o,M[j]!.output[i]+d);
            d := d+Size(AlphabetOfFRObject(M[j]));
        od;
        Add(trans,t);
        Add(out,o);
    od;
    x := FRMachineNC(FRMFamily([1..d]),N!.free,trans,out);
    SET_NAME@(M,"(+)",x);
    return x;
end);

InstallMethod(TensorProductOp, "(FR) for two FR machines",
        [IsList, IsFRMachine and IsFRMachineStdRep],
        function(M, N)
    local trans, out, t, o, i, j, x, a, b, alphabet, s;

    while ForAny(M,x->x!.free<>N!.free) do
        Error("All machines should have same stateset");
    od;

    alphabet := Cartesian(List(M,AlphabetOfFRObject));

    trans := [];
    out := [];
    for i in [1..Length(N!.output)] do
        t := [];
        o := [];
        for a in alphabet do
            b := [];
            s := i;
            for j in [1..Length(M)] do
                Add(b,Output(M[j],s,a[j]));
                s := Transition(M[j],s,a[j]);
            od;
            Add(o,Position(alphabet,b));
            Add(t,s);
        od;
        Add(trans,t);
        Add(out,o);
    od;
    x := FRMachineNC(FRMFamily([1..Length(alphabet)]),N!.free,trans,out);
    SET_NAME@(M,"(*)",x);
    return x;
end);

InstallMethod(DirectSumOp, "(FR) for two FR machines",
        [IsList, IsFRMachine and IsFRMachineStdRep],
        function(M, N)
    local c, gen, trans, out, t, o, i, j, d, alph, shift, sum;

    c := CallFuncList(LARGESTDENOMINATOR@,M);
    gen := Remove(c);

    d := 0;
    alph := []; shift := [];
    for i in [1..Length(M)] do
        Add(alph,AlphabetOfFRObject(M[i]));
        Add(shift, [d+1..d+Length(alph[i])]);
        d := d+Length(alph[i]);
    od;

    trans := [];
    out := [];
    for i in [1..Length(M)] do
        for j in gen[i] do
            t := ListWithIdenticalEntries(d,j^c[i]);
            t{shift[i]} := List(alph[i],a->Transition(M[i],j,a)^c[i]);
            o := [1..d];
            o{shift[i]} := shift[i]{Output(M[i],j)};
            Add(trans,t);
            Add(out,o);
        od;
    od;
    sum := FRMachineNC(FRMFamily([1..d]),Range(c[1]),trans,out);
    SetCorrespondence(sum,c);
    SET_NAME@(M,"#",sum);
    return sum;
end);

InstallMethod(DirectProductOp, "(FR) for two FR machines",
        [IsList, IsFRMachine and IsFRMachineStdRep],
        function(M, N)
    local c, gen, trans, out, t, o, i, j, a, b, product, alphabet;

    c := CallFuncList(LARGESTDENOMINATOR@,M);
    gen := Remove(c);

    alphabet := Cartesian(List(M,AlphabetOfFRObject));

    trans := [];
    out := [];
    for i in [1..Length(M)] do
        for j in gen[i] do
            t := [];
            o := [];
            for a in alphabet do
                b := ShallowCopy(a);
                b[i] := Output(M[i],j,a[i]);
                Add(o,Position(alphabet,b));
                Add(t,Transition(M[i],j,a[i])^c[i]);
            od;
            Add(trans,t);
            Add(out,o);
        od;
    od;
    product := FRMachineNC(FRMFamily([1..Length(alphabet)]),Range(c[1]),trans,out);
    SetCorrespondence(product,c);
    SET_NAME@(M,"x",product);
    return product;
end);

InstallMethod(TreeWreathProduct, "for two FR machines",
        [IsFRMachine and IsFRMachineStdRep, IsFRMachine and IsFRMachineStdRep,
         IsObject, IsObject],
        function(g,h,x0,y0)
    local c, gen, m, trans, out, t, o, one, i, j, alphabet;

    alphabet := Cartesian(AlphabetOfFRObject(g),AlphabetOfFRObject(h));
    while not [x0,y0] in alphabet do
        Error("(x0,y0) must be in the product of the machines' alphabets");
    od;
    c := LARGESTDENOMINATOR@(g,h,g,Zero(g));
    gen := Remove(c);
    if gen[4]=[] then
        one := One(Range(c[1]));
    else
        one := gen[4][1]^c[4];
    fi;

    trans := [];
    out := [];
    for i in [1..Length(gen[1])] do
        t := [];
        o := [];
        for j in alphabet do
            if j=[x0,y0] then
                Add(t,gen[1][i]^c[1]);
            elif j[2]=y0 then
                Add(t,gen[3][i]^c[3]);
            else
                Add(t,one);
            fi;
            Add(o,Position(alphabet,j));
        od;
        Add(trans,t);
        Add(out,o);
    od;
    for i in [1..Length(gen[2])] do
        t := [];
        o := [];
        for j in alphabet do
            if j[1]=x0 then
                Add(t,Transition(h,gen[2][i],j[2])^c[2]);
            else
                Add(t,one);
            fi;
            Add(o,Position(alphabet,[j[1],Output(h,gen[2][i],j[2])]));
        od;
        Add(trans,t);
        Add(out,o);
    od;
    for i in [1..Length(gen[3])] do
        t := [];
        o := [];
        for j in alphabet do
            if j[2]=y0 then
                Add(t,Transition(g,gen[3][i],j[1])^c[3]);
                Add(o,Position(alphabet,[Output(g,gen[3][i],j[1]),y0]));
            else
                Add(t,one);
                Add(o,Position(alphabet,j));
            fi;
        od;
        Add(trans,t);
        Add(out,o);
    od;
    Add(trans,ListWithIdenticalEntries(Length(alphabet),one));
    Add(out,[1..Length(alphabet)]);

    m := Minimized(FRMachineNC(FRMFamily([1..Length(alphabet)]),Range(c[1]),trans,out));
    m!.Correspondence := List(c{[1..2]},x->x*Correspondence(m));
    SET_NAME@([g,h],"~",m);
    return m;
end);
#############################################################################

#############################################################################
##
#M \=(FRMachine, FRMachine)
##
InstallMethod(\=, "(FR) for two FR machines",
        IsIdenticalObj,
        [IsFRMachine and IsFRMachineStdRep, IsFRMachine and IsFRMachineStdRep],
        function(left, right)
    local i, j;
    if left!.output <> right!.output then
        return false;
    elif Length(left!.transitions) <> Length(right!.transitions) then
        return false;
    fi;
    for i in [1..Length(left!.transitions)] do
        for j in AlphabetOfFRObject(left) do
            if LetterRepAssocWord(left!.transitions[i][j])<>LetterRepAssocWord(right!.transitions[i][j]) then return false; fi;
        od;
    od;
    return true;
end);
#############################################################################

#############################################################################
##
#M \<(FRMachine, FRMachine)
##
InstallMethod(\<, "(FR) for two FR machines",
        IsIdenticalObj,
        [IsFRMachine and IsFRMachineStdRep, IsFRMachine and IsFRMachineStdRep],
        function(left, right)
    local i, j, wl, wr;
    if left!.output <> right!.output then
        return left!.output < right!.output;
    elif Length(left!.transitions) <> Length(right!.transitions) then
        return Length(left!.transitions) < Length(right!.transitions);
    fi;
    for i in [1..Length(left!.transitions)] do
        for j in AlphabetOfFRObject(left) do
            wl := LetterRepAssocWord(left!.transitions[i][j]);
            wr := LetterRepAssocWord(right!.transitions[i][j]);
            if wl<>wr then return wl<wr; fi;
        od;
    od;
    return false; # they're equal
end);
#############################################################################

#############################################################################
##
#A WreathRecursion(FRMachine)
#O Output(Machine, State)
#O Transition(Machine, State, Letter)
##
InstallMethod(Output, "(FR) for an FR machine",
        [IsGroupFRMachine and IsFRMachineStdRep],
        function(M)
    local image;
    image := List(M!.output,PermList);
    return GroupHomomorphismByImages(StateSet(M),Group(image),image);
end);

InstallMethod(Output, "(FR) for an FR machine",
        [IsMonoidFRMachine and IsFRMachineStdRep],
        function(M)
    local image;
    image := List(M!.output,TransformationList);
    return SemigroupHomomorphismByImagesNC(StateSet(M),Monoid(image),image);
end);

InstallMethod(Output, "(FR) for an FR machine",
        [IsSemigroupFRMachine and IsFRMachineStdRep],
        function(M)
    local image;
    image := List(M!.output,TransformationList);
    return SemigroupHomomorphismByImagesNC(StateSet(M),Semigroup(image),image);
end);

InstallMethod(Output, "(FR) for an FR machine and a state expressed as an integer",
        [IsFRMachine and IsFRMachineStdRep, IsInt],
        function(M, i)
    if i > 0 then
        return M!.output[i];
    elif i = 0 then
        return AlphabetOfFRObject(M);
    else
        return INVERSE@(M!.output[-i]);
    fi;
end);

InstallMethod(Output, "(FR) for an FR machine and a state expressed as a word",
        [IsFRMachine and IsFRMachineStdRep, IsAssocWord],
        function(M, w)
    local perm, i;
    perm := AlphabetOfFRObject(M);
    for i in LetterRepAssocWord(w) do
        if i > 0 then
            perm := M!.output[i]{perm};
        else
            perm := INVERSE@(M!.output[-i]){perm};
        fi;
    od;
    return perm;
end);

InstallMethod(Output, "(FR) for an FR machine and a state expressed as a list",
        [IsFRMachine, IsList],
        function(M, l)
    local perm, i;
    perm := AlphabetOfFRObject(M);
    for i in l do
        perm := Output(M,i){perm};
    od;
    return perm;
end);

InstallMethod(Output, "(FR) for an FR machine, a state and a letter",
        [IsFRMachine, IsObject, IsObject],
        function(M, s, a)
    return Output(M,s)[a];
end);

InstallMethod(Transition, "(FR) for an FR machine, a state expressed as an integer, and an input",
        [IsFRMachine and IsFRMachineStdRep, IsInt, IsPosInt],
        function(M, i, p)
    if i > 0 then
        return M!.transitions[i][p];
    else
        return M!.transitions[-i][PREIMAGE@(M!.output[-i],p)];
    fi;
end);

InstallMethod(Transitions, "(FR) for an FR machine and a state expressed as an integer",
        [IsFRMachine and IsFRMachineStdRep, IsInt],
        function(M, i)
    if i > 0 then
        return M!.transitions[i];
    else
        return M!.transitions[-i]{INVERSE@(M!.output[-i])};
    fi;
end);

BindGlobal("FRMTRANSITION@", function(M,l,p)
    local w, i;
    if IsMonoid(M!.free) then
        w := One(M!.free);
    else
        w := fail;
    fi;
    for i in l do
        if i > 0 then
            if w=fail then
                w := M!.transitions[i][p];
            else
                w := w*M!.transitions[i][p];
            fi;
            p := M!.output[i][p];
        else
            p := PREIMAGE@(M!.output[-i],p);
            if w=fail then
                w := M!.transitions[-i][p]^-1;
            else
                w := w/M!.transitions[-i][p];
            fi;
        fi;
    od;
    return w;
end);

InstallMethod(Transition, "(FR) for an FR machine, a state expressed as a list, and an input",
        [IsFRMachine and IsFRMachineStdRep, IsList, IsPosInt],
        FRMTRANSITION@);

InstallMethod(Transition, "(FR) for an FR machine, a state expressed as a word, and an input",
        [IsFRMachine and IsFRMachineStdRep, IsAssocWord, IsPosInt],
        function(M, v, p)
    return FRMTRANSITION@(M,LetterRepAssocWord(v),p);
end);

InstallMethod(Transition, "(FR) for an FR machine, a state, and a list of letters",
        [IsFRMachine, IsObject, IsList],
        function(M, s, l)
    local i, t;
    t := s;
    for i in l do t := Transition(M, t, i); od;
    return t;
end);

InstallMethod(Transitions, "(FR) for an FR machine and a state expressed as a list",
        [IsFRMachine and IsFRMachineStdRep, IsList],
        function(M,w)
    return WreathRecursion(M)(w)[1];
end);

InstallMethod(Transitions, "(FR) for an FR machine, a state expressed as a word, and an input",
        [IsFRMachine and IsFRMachineStdRep, IsAssocWord],
        function(M, w)
    return WreathRecursion(M)(w)[1];
end);

InstallMethod(WreathRecursion, "(FR) for an FR machine",
        [IsFRMachine and IsFRMachineStdRep],
        function(M)
    local output, alphabet, ones, transitions, onetrans;
    output := M!.output;
    alphabet := AlphabetOfFRObject(M);
    if IsMonoid(M!.free) then
        ones := List(alphabet,x->One(M!.free));
    else
        ones := List(alphabet,x->fail);
    fi;
    onetrans := AlphabetOfFRObject(M);
    transitions := M!.transitions;
    return function(w)
        local vector, perm, i, j;
        vector := ShallowCopy(ones);
        perm := onetrans;
        if IsAssocWord(w) then w := LetterRepAssocWord(w); fi;
        for i in w do
            if i > 0 then
                if vector[1]=fail then
                    vector := ShallowCopy(transitions[i]);
                else
                    for j in alphabet do
                        vector[j] := vector[j]*transitions[i][perm[j]];
                    od;
                fi;
                perm := output[i]{perm};
            else
                perm := INVERSE@(output[-i]){perm};
                if vector[1]=fail then
                    vector := List(transitions[-i]{perm},Inverse);
                else
                    for j in alphabet do
                        vector[j] := vector[j]/transitions[-i][perm[j]];
                    od;
                fi;
            fi;
        od;
        return [vector,perm];
    end;
end);

InstallMethod(VirtualEndomorphism, "(FR) for a group FR machine and a vertex",
        [IsGroupFRMachine,IsObject],
        function(M,v)
    local G, H;
    G := StateSet(M);
    H := Stabilizer(G,v,function(w,g) return w^FRElement(M,g); end);
    return GroupHomomorphismByImages(H,G,GeneratorsOfGroup(H),List(GeneratorsOfGroup(H),x->Transition(M,x,v)));
end);
#############################################################################

#############################################################################
##
#A FRMachineRWS
##
InstallMethod(FRMachineRWS, "(FR) for an FR machine",
        [IsFRMachine and IsFRMachineStdRep],
        function(M)
    local output, alphabet, transitions, inverse, iso, fpmonoid, gens, mgens, rws;

    if IsGroupFRMachine(M) then
        iso := IsomorphismFpMonoidInversesFirst(M!.free);
        mgens := List(GeneratorsOfMonoid(Range(iso)),x->PreImage(iso,x));
        rws := rec(rws := KnuthBendixRewritingSystem(Range(iso)),
                   letterrep := w->LetterRepAssocWord(UnderlyingElement(w^iso)),
                   letterunrep := w->Product(mgens{w},One(M!.free)));
        gens := List(GeneratorsOfMonoid(Range(iso)),
                     w->FRElement(M,PreImagesRepresentativeNC(iso,w)));
        inverse := List(gens,w->rws.letterrep(InitialState(w)^-1)[1]);
        rws.cyclicallyreduce := function(w)
            local i, j;
            i := 1; j := Length(w);
            while i<j and w[i]=inverse[w[j]] do
                i := i+1; j := j-1;
            od;
            if i=1 then return w; else return w{[i..j]}; fi;
        end;
    else
        rws := rec(rws := KnuthBendixRewritingSystem(M!.free/[]),
                   letterrep := LetterRepAssocWord);
        if IsMonoidFRMachine(M) then
            rws.letterunrep := w->Product(GeneratorsOfMonoid(M!.free){w},One(M!.free));
        else
            rws.letterunrep := w->Product(GeneratorsOfSemigroup(M!.free){w});
        fi;
        gens := List(GeneratorsOfFRMachine(M),w->FRElement(M,w));
    fi;
    output := List(gens,Output);
    alphabet := AlphabetOfFRObject(M);
    transitions := List(gens,w->List(alphabet,
                           a->rws.letterrep(Transition(w,a))));
    if ValueOption("fr_maxlen")<>fail then
        rws.maxlen := ValueOption("fr_maxlen");
    else
        rws.maxlen := 5; # do not add rules longer than that -- too slow
    fi;
    rws.modified := true; # whether the true rules rws.tzrules and
    # the temporary rules rws.rws!.tzrules are in sync
    rws.pi := function(w)
        local vector, perm, i, j;
        vector := List(alphabet,x->[]);
        perm := alphabet;
        for i in w do
            for j in alphabet do
                Append(vector[j],transitions[i][perm[j]]);
            od;
            perm := output[i]{perm};
        od;
        return [vector,perm];
    end;
    rws.reduce := w->ReduceLetterRepWordsRewSys(rws.rws!.tzrules,w);
    rws.addsgrule := function(l,r,short)
        local ll, lr;
        ll := Length(l); lr := Length(r);
        if short and (ll>rws.maxlen or lr>rws.maxlen) then return; fi;
        rws.modified := true;
        if ll>lr or (ll=lr and l>r) then
            Info(InfoFR,3,"# Added rule ",l," -> ",r);
            AddRuleReduced(rws.rws,[ShallowCopy(l),ShallowCopy(r)]);
        else
            Info(InfoFR,3,"# Added rule ",r," -> ",l);
            AddRuleReduced(rws.rws,[ShallowCopy(r),ShallowCopy(l)]);
        fi;
    end;
    if IsGroupFRMachine(M) then
        rws.addgprule := function(w,short)
            local i, l, ll, left, right;
            l := Length(w);
            if short and l>2*rws.maxlen then return; fi;
            Info(InfoFR,3,"# Adding group rule ",w);
            rws.addsgrule(w,[],short);
            rws.addsgrule(inverse{w{[l,l-1..1]}},[],short);
            ll := QuoInt(l,2);
            left := w{[1..ll]};
            right := inverse{w{[l,l-1..ll+1]}};
            for i in [1..l] do
                rws.addsgrule(left,right,short);
                Add(left,inverse[Remove(right)]);
                if IsOddInt(l) then
                    rws.addsgrule(left,right,short);
                fi;
                Add(right,inverse[Remove(left,1)],1);
            od;
        end;
    fi;
    rws.commit := function()
        if rws.modified then
            Info(InfoFR,3,"# Committed rules ",rws.rws!.tzrules);
            rws.tzrules := StructuralCopy(rws.rws!.tzrules);
            rws.modified := false;
        fi;
    end;
    rws.restart := function()
        if rws.modified then
            Info(InfoFR,3,"# Restarting with fresh rules ",rws.tzrules);
            rws.rws!.tzrules := StructuralCopy(rws.tzrules);
            rws.modified := false;
        fi;
    end;
    rws.commit();
    return rws;
end);

InstallGlobalFunction(NewFRMachineRWS, # "(FR) will restart with fresh rules",
        function(M)
    local rws;
    rws := FRMachineRWS(M);
    rws.restart();
    return rws;
end);
#############################################################################

############################################################################
##
#O States(FRMachine)
##
InstallMethod(States, "(FR) for an FR machine and an element",
        [IsFRMachine, IsMultiplicativeElement],
        function(M,x)
    return States(M,[x]);
end);

InstallOtherMethod(States, "(FR) for an FR machine and a list of elements",
        [IsFRMachine, IsMultiplicativeElementCollection],
        function(M,L)
    local states, i, x, stateset;
    states := ShallowCopy(L);
    stateset := Set(states);
    i := 1;
    while i <= Length(states) do
        for x in Transitions(M,states[i]) do
            if not x in stateset then
                Add(states,x);
                AddSet(stateset,x);
            fi;
        od;
        i := i+1;
        if RemInt(i,100)=0 then
            Info(InfoFR, 2, "The states contain at least ", states);
        fi;
    od;
    return states;
end);

InstallMethod(FixedStates, "(FR) for an FR machine and an element",
        [IsFRMachine, IsMultiplicativeElement],
        function(M,x)
    return FixedStates(M,[x]);
end);

InstallMethod(FixedStates, "(FR) for a list of FR elements",
        [IsFRMachine, IsMultiplicativeElementCollection],
        function(M,L)
    local states, alphabet, i, x, addstates, stateset;
    states := [];
    stateset := [];
    alphabet := AlphabetOfFRObject(M);
    addstates := function(x)
        local i, o, t;
        o := Output(M,x);
        t := Transitions(M,x);
        for i in [1..Length(alphabet)] do
            if o[i]=alphabet[i] and not t[i] in stateset then
                Add(states,t[i]);
                AddSet(stateset,t[i]);
            fi;
        od;
    end;
    for x in L do addstates(x); od;
    i := 1;
    while i <= Length(states) do
        addstates(states[i]);
        i := i+1;
        if RemInt(i,100)=0 then
            Info(InfoFR, 2, "The fixed states contain at least ", states);
        fi;
    od;
    return states;
end);

InstallMethod(LimitStates, "(FR) for an FR machine and element",
        [IsFRMachine,IsMultiplicativeElement],
        function(M,x)
    return LimitStates(M,[x]);
end);

InstallMethod(LimitStates, "(FR) for an FR machine and a list of elements",
        [IsFRMachine,IsMultiplicativeElementCollection],
        function(M,L)
    local s, d, S, oldS;
    s := Set(States(M,L));
    d := List(s,w->BlistList([1..Length(s)],List(DecompositionOfFRElement(w)[1],x->Position(s,x))));
    S := BlistList([1..Length(s)],[1..Length(s)]);
    repeat
        oldS := S;
        S := UnionBlist(ListBlist(d,S));
    until oldS=S;
    return ListBlist(s,S);
end);

InstallMethod(CoverNucleus, "(FR) for an FR machine",
        [IsFRMachine],
        function(M)
    local s, news, olds, gens, g, h;

    gens := Set(GeneratorsOfFRMachine(M));
    news := gens;
    s := [];

    while true do
        olds := ShallowCopy(s);
        UniteSet(s,LimitStates(M,news));
        if Length(s)=Length(olds) then
            return s;
        fi;

        news := [];
        for g in Difference(s,olds) do
            if g in FixedStates(g) and Order(g)=infinity then
                return fail;
            fi;
            for h in gens do AddSet(news,g*h); od;
        od;
        Info(InfoFR, 2, "Nucleus: The nucleus contains at least ",s);
    od;
end);

#############################################################################
##
#M StructuralGroup(FRMachine)
#M StructuralSemigroup(FRMachine)
#M StructuralMonoid(FRMachine)
##
InstallMethod(StructuralGroup, "(FR) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    local ggens, wgens, f, fggens, fwgens, phi;
    ggens := GeneratorsOfFRMachine(M);
    wgens := [1..Size(AlphabetOfFRObject(M))];
    f := FreeGroup(Concatenation(List(ggens,String),List(wgens,String)));
    fggens := GeneratorsOfGroup(f){[1..Length(ggens)]};
    fwgens := GeneratorsOfGroup(f){[Length(ggens)+1..Length(ggens)+Length(wgens)]};
    phi := GroupHomomorphismByImagesNC(StateSet(M),f,ggens,fggens);
    return f / List(Cartesian([1..Length(ggens)],wgens),
                 p->fggens[p[1]]*fwgens[Output(M,p[1],p[2])]/Transition(M,p[1],p[2])^phi/fwgens[p[2]]);
end);

InstallMethod(StructuralMonoid, "(FR) for a monoid FR machine",
        [IsMonoidFRMachine],
        function(M)
    local ggens, wgens, f, fggens, fwgens;
    ggens := GeneratorsOfFRMachine(M);
    wgens := [1..Size(AlphabetOfFRObject(M))];
    f := FreeMonoid(Concatenation(List(ggens,String),List(wgens,String)));
    fggens := GeneratorsOfMonoid(f){[1..Length(ggens)]};
    fwgens := GeneratorsOfMonoid(f){[Length(ggens)+1..Length(ggens)+Length(wgens)]};
    return f / List(Cartesian([1..Length(ggens)],wgens),
                 p->[fggens[p[1]]*fwgens[Output(M,p[1],p[2])],fwgens[p[2]]*MappedWord(Transition(M,p[1],p[2]),ggens,fggens)]);
end);

InstallMethod(StructuralSemigroup, "(FR) for a semigroup FR machine",
        [IsSemigroupFRMachine],
        function(M)
    local ggens, wgens, f, fggens, fwgens;
    ggens := GeneratorsOfFRMachine(M);
    wgens := [1..Size(AlphabetOfFRObject(M))];
    f := FreeSemigroup(Concatenation(List(ggens,String),List(wgens,String)));
    fggens := GeneratorsOfSemigroup(f){[1..Length(ggens)]};
    fwgens := GeneratorsOfSemigroup(f){[Length(ggens)+1..Length(ggens)+Length(wgens)]};
    return f / List(Cartesian([1..Length(ggens)],wgens),
                 p->[fggens[p[1]]*fwgens[Output(M,p[1],p[2])],fwgens[p[2]]*MappedWord(Transition(M,p[1],p[2]),ggens,fggens)]);
end);
#############################################################################

#############################################################################
##
#M AsGroupFRMachine
#M AsMonoidFRMachine
#M AsSemigroupFRMachine
##
InstallMethod(AsGroupFRMachine, "(FR) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    SetCorrespondence(M,IdentityMapping(M!.free));
    return M;
end);

BindGlobal("ASGROUPFRMACHINE@", function(M)
    local f, N, h, s;
    if not ForAll(M!.output,ISINVERTIBLE@) then return fail; fi;
    s := GeneratorsOfFRMachine(M);
    f := FreeGroup(Length(s));
    h := MagmaHomomorphismByImagesNC(M!.free,f,GeneratorsOfGroup(f));
    N := FRMachineNC(FamilyObj(M),f,List(M!.transitions,r->List(r,w->w^h)),M!.output);
    SetCorrespondence(N,h);
    return N;
end);
InstallMethod(AsGroupFRMachine, "(FR) for a monoid FR machine",
        [IsMonoidFRMachine],
        ASGROUPFRMACHINE@);
InstallMethod(AsGroupFRMachine, "(FR) for a semigroup FR machine",
        [IsSemigroupFRMachine],
        ASGROUPFRMACHINE@);

InstallMethod(AsMonoidFRMachine, "(FR) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    local f, N, h, sM, sN, sNinv, trans, out, o, i;
    sM := GeneratorsOfGroup(M!.free);
    f := FreeMonoid(2*Length(sM));
    sN := GeneratorsOfMonoid(f){[1..Length(sM)]};
    sNinv := GeneratorsOfMonoid(f){[Length(sM)+1..2*Length(sM)]};
    h := MagmaHomomorphismByFunctionNC(M!.free,f,function(w)
        local r, i;
        r := [];
        for i in LetterRepAssocWord(w) do
            if i>0 then Add(r,i); else Add(r,Length(sM)-i); fi;
        od;
        return AssocWordByLetterRep(FamilyObj(sN[1]),r);
    end);
    trans := List(M!.transitions,r->List(r,w->w^h));
    out := ShallowCopy(M!.output);
    for i in [1..Length(M!.transitions)] do
        o := INVERSE@(M!.output[i]);
        Add(trans,List(M!.transitions[i],w->(w^-1)^h){o});
        Add(out,o);
    od;
    N := FRMachineNC(FamilyObj(M),f,trans,out);
    SetCorrespondence(N,h);
    return N;
end);

InstallMethod(AsMonoidFRMachine, "(FR) for a monoid FR machine",
        [IsMonoidFRMachine],
        function(M)
    SetCorrespondence(M,IdentityMapping(M!.free));
    return M;
end);

InstallMethod(AsMonoidFRMachine, "(FR) for a semigroup FR machine",
        [IsSemigroupFRMachine],
        function(M)
    local f, N, h, s;
    s := GeneratorsOfSemigroup(M!.free);
    f := FreeMonoid(Length(s));
    h := MagmaHomomorphismByImagesNC(M!.free,f,GeneratorsOfMonoid(f));
    N := FRMachineNC(FamilyObj(M),f,List(M!.transitions,r->List(r,w->w^h)),M!.output);
    SetCorrespondence(N,h);
    return N;
end);

InstallMethod(AsSemigroupFRMachine, "(FR) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    local f, N, h, sM, sN, sNinv, one, trans, out, o, i;
    sM := GeneratorsOfGroup(M!.free);
    f := FreeSemigroup(2*Length(sM)+1);
    sN := GeneratorsOfSemigroup(f){[1..Length(sM)]};
    sNinv := GeneratorsOfSemigroup(f){[Length(sM)+1..2*Length(sM)]};
    one := GeneratorsOfSemigroup(f)[2*Length(sM)+1];
    h := MagmaHomomorphismByFunctionNC(M!.free,f,function(w)
        local r, i;
        r := [];
        if IsOne(w) then
            return one;
        else
            for i in LetterRepAssocWord(w) do
                if i>0 then Add(r,i); else Add(r,Length(sM)-i); fi;
            od;
            return AssocWordByLetterRep(FamilyObj(one),r);
        fi;
    end);
    trans := List(M!.transitions,r->List(r,w->w^h));
    out := ShallowCopy(M!.output);
    for i in [1..Length(M!.transitions)] do
        o := INVERSE@(M!.output[i]);
        Add(trans,List(M!.transitions[i],w->(w^-1)^h){o});
        Add(out,o);
    od;
    Add(trans,List(AlphabetOfFRObject(M),a->one)); # add an identity state
    Add(out,AlphabetOfFRObject(M));
    N := FRMachineNC(FamilyObj(M),f,trans,out);
    SetCorrespondence(N,h);
    return N;
end);

InstallMethod(AsSemigroupFRMachine, "(FR) for a monoid FR machine",
        [IsMonoidFRMachine],
        function(M)
    local f, N, h, sM, sN, one, trans, out, i;
    sM := GeneratorsOfMonoid(M!.free);
    f := FreeSemigroup(Length(sM)+1);
    sN := GeneratorsOfSemigroup(f){[1..Length(sM)]};
    one := GeneratorsOfSemigroup(f)[Length(sM)+1];
    h := MagmaHomomorphismByFunctionNC(M!.free,f,w->MAPPEDWORD@(w,sN,one));
    trans := List(M!.transitions,r->List(r,w->w^h));
    out := ShallowCopy(M!.output);
    Add(trans,List(AlphabetOfFRObject(M),a->one));
    Add(out,AlphabetOfFRObject(M));
    N := FRMachineNC(FamilyObj(M),f,trans,out);
    SetCorrespondence(N,h);
    return N;
end);

InstallMethod(AsSemigroupFRMachine, "(FR) for a semigroup FR machine",
        [IsSemigroupFRMachine],
        function(M)
    SetCorrespondence(M,IdentityMapping(M!.free));
    return M;
end);

BindGlobal("HOM2MACHINE@", function(f,tester,g)
    local s;
    s := Source(f);
    if not tester(s) or s<>Range(f) then
        return fail;
    fi;
    return FRMachineNC(FRMFamily([1]),s,List(g(s),x->[x^f]),List(g(s),x->[1]));
end);
    
InstallMethod(AsGroupFRMachine, "(FR) for a group homomorphism",
        [IsGroupHomomorphism],
        f->HOM2MACHINE@(f,IsFreeGroup,GeneratorsOfGroup));

InstallMethod(AsMonoidFRMachine, "(FR) for a monoid homomorphism",
        [IsMagmaHomomorphism],
        f->HOM2MACHINE@(f,IsFreeMonoid,GeneratorsOfMonoid));

InstallMethod(AsSemigroupFRMachine, "(FR) for a semigroup homomorphism",
        [IsMagmaHomomorphism],
        f->HOM2MACHINE@(f,IsFreeSemigroup,GeneratorsOfSemigroup));
#############################################################################

#############################################################################
##
#M Minimized(FRMachine)
##
BindGlobal("MINIMIZERWS_MAKERULES@", function(rws,p)
    # p is a tuple [generators,inverses,isone?,rules]
    # this command recomputes the rules
    local i, l;
    p[4] := [];
    if p[3] then
        for i in p[1] do
            Add(p[4],[[i],[]]);
        od;
    else
        l := p[1][Length(p[1])];
        for i in [1..Length(p[1])-1] do
            Add(p[4],[[p[1][i]],[l]]);
        od;
        if l in p[2] then
            Add(p[4],[[l,l],[]]);
        fi;
    fi;
end);

BindGlobal("MINIMIZERWS@", function(M)
    local rws, gens, h, i, j, si, p, part, newpart, changed;

    rws := NewFRMachineRWS(M);
    if IsBound(rws.partition) then
        return rws;
    fi;
    if IsSemigroupFRMachine(M) then
        gens := GeneratorsOfSemigroup(M!.free);
    else
        gens := GeneratorsOfMonoid(M!.free);
    fi;
    gens := Filtered(gens,x->rws.letterrep(x)=rws.reduce(rws.letterrep(x)));
    i := List(gens,x->Output(M,x));
    si := Set(i);
    part := List(si,x->[[],[],HasIsBuiltFromMonoid(rws.rws) and IsBuiltFromMonoid(rws.rws) and ISONE@(x)]);
    for j in [1..Length(i)] do
        p := Position(si,i[j]);
        Add(part[p][1],rws.letterrep(gens[j])[1]);
        if IsGroupFRMachine(M) then
            Add(part[p][2],rws.letterrep(gens[j]^-1)[1]);
        else
            Add(part[p][2],fail);
        fi;
    od;
    for p in part do
        SortParallel(p[1],p[2]);
        MINIMIZERWS_MAKERULES@(rws,p);
    od;

    changed := true;
    while changed do
        #Info(InfoFR,1,"New parts: ",part);
        changed := false;
        rws.rws!.tzrules := Concatenation(rws.tzrules,Concatenation(List(part,p->p[4])));
        for h in [1..Length(part)] do
            i := List(part[h][1],x->List(rws.pi([x])[1],rws.reduce));
            si := Set(i);
            if Length(si)>1 then
                changed := true;
                newpart := List(si,x->[[],[],part[h][3] and ForAll(x,IsEmpty)]);
                for j in [1..Length(i)] do
                    p := Position(si,i[j]);
                    Add(newpart[p][1],part[h][1][j]);
                    Add(newpart[p][2],part[h][2][j]);
                od;
                for p in newpart do
                    MINIMIZERWS_MAKERULES@(rws,p);
                od;
                Append(part,newpart);
                part[h] := Remove(part);
                break;
            elif part[h][3] and not ForAll(si[1],IsEmpty) then
                changed := true;
                part[h][3] := false;
                MINIMIZERWS_MAKERULES@(rws,part[h]);
                break;
            fi;
        od;
    od;
    for p in part do
        p[3] := p[3] and ForAll(p[1],x->ForAll(rws.pi([x])[1],x->rws.reduce(x)=[]));
        MINIMIZERWS_MAKERULES@(rws,p);
    od;
    rws.rws!.tzrules := Concatenation(rws.tzrules,Concatenation(List(part,p->p[4])));
    rws.modified := true;
    rws.commit();
    rws.partition := part;
    return rws;
end);

InstallMethod(Minimized, "(FR) for a group/monoid/semigroup FR machine",
        [IsFRMachine and IsFRMachineStdRep],
        function(M)
    local rws, gens, gensimg, i, ri, red, free, freegens, one, out, trans, map;
    rws := MINIMIZERWS@(M);
    gens := GeneratorsOfFRMachine(M);
    i := List(gens,rws.letterrep);
    red := Filtered(i,x->rws.reduce(x)=x);
    if i=red then
        M := COPYFRMACHINE@(M);
        SetCorrespondence(M,IdentityMapping(M!.free));
        return M;
    fi;
    if IsGroupFRMachine(M) then
        free := FreeGroup(Length(red));
        freegens := GeneratorsOfGroup(free);
        one := One(free);
    elif IsMonoidFRMachine(M) then
        free := FreeMonoid(Length(red));
        freegens := GeneratorsOfMonoid(free);
        one := One(free);
    elif IsSemigroupFRMachine(M) then
        free := FreeSemigroup(Length(red));
        freegens := GeneratorsOfSemigroup(free);
        one := fail;
    fi;
    gensimg := [];
    for i in gens do
        ri := rws.reduce(rws.letterrep(i));
        if ri=[] then
            Add(gensimg,One(free));
        elif ri in red then
            Add(gensimg,freegens[Position(red,ri)]);
        else
            Add(gensimg,freegens[Position(red,rws.reduce(rws.letterrep(i^-1)))]^-1);
        fi;
    od;
    map := MagmaHomomorphismByImagesNC(M!.free,free,gensimg);
    i := List(red,rws.pi);

    out := List(i,p->p[2]);
    trans := [];
    for i in i do
        Add(trans,List(i[1],w->rws.letterunrep(rws.reduce(w))^map));
    od;
    i := FRMachineNC(FamilyObj(M),free,trans,out);
    SetCorrespondence(i,map);
    return i;
end);
#############################################################################

#############################################################################
##
#M SubFRMachine(FRMachine,FRMachine)
##
InstallMethod(SubFRMachine, "(FR) for two group/monoid/semigroup FR machines",
        [IsFRMachine and IsFRMachineStdRep,IsFRMachine and IsFRMachineStdRep],
        function(M,N)
    local rws, S, Mgens, Ngens, Mletter, Nred;

    if AlphabetOfFRObject(M)<>AlphabetOfFRObject(N) then
        return fail;
    elif IsIdenticalObj(M,N) then
        return IdentityMapping(M!.free);
    elif M=N then
        if IsGroupFRMachine(M) then
            return GroupHomomorphismByImages(N!.free,M!.free,GeneratorsOfGroup(N!.free),GeneratorsOfGroup(M!.free));
        else
            return MagmaHomomorphismByImagesNC(N!.free,M!.free,GeneratorsOfFRMachine(M));
        fi;
    fi;
    if (IsGroupFRMachine(N) and not IsGroupFRMachine(M)) or (IsMonoidFRMachine(N) and IsSemigroupFRMachine(M)) then
        return fail;
    fi;
    S := FRMMINSUM@(N,M);
    rws := MINIMIZERWS@(S);
    Mgens := GeneratorsOfSemigroup(M!.free);
    Ngens := GeneratorsOfFRMachine(N);
    Mletter := List(Mgens,x->rws.letterrep(x^Correspondence(S)[2]));
    Nred := List(Ngens,x->rws.reduce(rws.letterrep(x^Correspondence(S)[1])));
    if IsSubset(Mletter,Nred) then
        Mgens := List(Nred,x->Mgens[Position(Mletter,x)]);
        return MagmaHomomorphismByImagesNC(N!.free,M!.free,Mgens);
    else
        return fail;
    fi;
end);

InstallMethod(SubFRMachine, "(FR) for a machine and a homomorphism",
        [IsFRMachine and IsFRMachineStdRep, IsMapping],
        function(M,f)
    local S, trans, out, i, pi, x;
    S := StateSet(M);
    while S<>Range(f) do
        Error("SubFRMachine: range and stateset must be the same\n");
    od;
    while not IsFreeGroup(Source(f)) do
        Error("SubFRMachine: source must be a free group\n");
    od;
    pi := WreathRecursion(M);
    trans := [];
    out := [];
    for i in GeneratorsOfGroup(Source(f)) do
        x := pi(i^f);
        x[1] := List(x[1],g->PreImagesRepresentativeNC(f,g));
        if fail in x[1] then return fail; fi;
        Add(trans,x[1]);
        Add(out,x[2]);
    od;
    x := FRMachineNC(FamilyObj(M),Source(f),trans,out);
    if HasAddingElement(M) then
        i := PreImagesRepresentativeNC(f,InitialState(AddingElement(M)));
        if i<>fail then
            SetAddingElement(x,FRElement(x,i));
        fi;
    fi;
    return x;
end);
#############################################################################

################################################################
# change basis of FR machine
BindGlobal("CHANGEFRMACHINEBASIS@", function(M,l,p)
    local trans, i, d, newM;
    d := Size(AlphabetOfFRObject(M));
    while Length(l)<>d or not ForAll(l,x->x in StateSet(M)) do
        Error("Invalid base change ",l,"\n");
    od;
    while LargestMovedPoint(p)>d do
	Error("Invalid permutation ",p,"\n");
    od;
    trans := [];
    for i in [1..Length(M!.transitions)] do
        Add(trans,Permuted(List(AlphabetOfFRObject(M),a->l[a]^-1*M!.transitions[i][a]*l[M!.output[i][a]]),p));
    od;
    newM := FRMachineNC(FamilyObj(M),StateSet(M),trans,List(M!.output,r->ListPerm(PermList(r)^p,d)));
    return newM;
end);

InstallMethod(ChangeFRMachineBasis, "(FR) for a group FR machine and a list",
        [IsGroupFRMachine, IsCollection],
        function(M,l)
    return ChangeFRMachineBasis(M,l,());
end);	
InstallMethod(ChangeFRMachineBasis, "(FR) for a group FR machine and a permutation",
        [IsGroupFRMachine, IsPerm],
        function(M,p)
    return ChangeFRMachineBasis(M,List(AlphabetOfFRObject(M),x->One(StateSet(M))),p);
end);	
InstallMethod(ChangeFRMachineBasis, "(FR) for a group FR machine, a list and a permutation",
        [IsGroupFRMachine, IsCollection, IsPerm],
    CHANGEFRMACHINEBASIS@);

InstallMethod(ChangeFRMachineBasis, "(FR) for an FR machine",
        [IsGroupFRMachine],
        function(M)
    local cycles, basis, s, t, u, v;
    
    # gather all permutation cycles
    cycles := [];
    for s in GeneratorsOfFRMachine(M) do
        for t in Cycles(PermList(Output(M,s)),AlphabetOfFRObject(M)) do
            if Length(t)>1 then
                Add(cycles,[s,t]);
            fi;
        od;
    od;
    
    basis := [];
    while cycles<>[] do
        # first cycle connected to the partial basis
        t := First([1..Length(cycles)],i->Number(cycles[i][2],i->IsBound(basis[i]))>0);
        if t=fail then
            # set up an anchor on the cycle
            basis[First(AlphabetOfFRObject(M),i->not IsBound(basis[i]))] := One(StateSet(M));
            continue;
        fi;
        t := Remove(cycles,t);
        # anchor on the cycle
        s := First(t[2],i->IsBound(basis[i]));
        u := s;
        repeat
            v := Output(M,t[1],u);
            if not IsBound(basis[v]) then
                basis[v] := LeftQuotient(Transition(M,t[1],u),basis[u]);
            fi;
            u := v;
        until u=s;
    od;
    return ChangeFRMachineBasis(M,basis,());
end);

BindGlobal("RIGHTACTMACHINE@", function(M,f)
    local S;
    S := StateSet(M);
    if S<>Source(f) or S<>Range(f) then
        Error("\*: source, range and stateset must be the same\n");
    fi;
    return FRMachineNC(FamilyObj(M),S,List(M!.transitions,r->List(r,x->x^f)),M!.output);
end);

InstallMethod(\*, "(FR) for an FR machine and a mapping",
        [IsFRMachine and IsFRMachineStdRep, IsMapping],
        RIGHTACTMACHINE@);

BindGlobal("LEFTACTMACHINE@", function(f,M)
    local S, trans, out, i, pi, x;
    S := StateSet(M);
    if S<>Source(f) or S<>Range(f) then
        Error("\*: source, range and stateset must be the same\n");
    fi;
    pi := WreathRecursion(M);
    trans := [];
    out := [];
    
    for i in [1..Length(M!.output)] do
        x := pi(GeneratorsOfFRMachine(M)[i]^f);
        Add(trans,x[1]);
        Add(out,x[2]);
    od;
    return FRMachineNC(FamilyObj(M),S,trans,out);
end);

InstallMethod(\*, "(FR) for a mapping and an FR machine",
        [IsMapping, IsFRMachine and IsFRMachineStdRep],
        LEFTACTMACHINE@);

BindGlobal("CONJACTMACHINE@", function(M,f)
    local S, newS, trans, out, i, pi, x, finv;
    S := StateSet(M);
    if S<>Source(f) then
        Error("\^: source and stateset must be the same\n");
    fi;
    newS := Range(f);
    pi := WreathRecursion(M);
    trans := [];
    out := [];
    finv := InverseGeneralMapping(f);
    if finv=fail then return fail; fi;
    for i in GeneratorsOfGroup(newS) do
        x := pi(ImagesRepresentative(finv,i));
        Add(trans,List(x[1],x->ImagesRepresentative(f,x)));
        Add(out,x[2]);
    od;
    return FRMachineNC(FamilyObj(M),newS,trans,out);
end);

InstallMethod(\^, "(FR) for a group FR machine and a mapping",
        [IsFRMachine and IsFRMachineStdRep, IsMapping],
        CONJACTMACHINE@);
################################################################

