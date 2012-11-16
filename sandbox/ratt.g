if not IsBound(VHStructure) then
DeclareAttribute("VHStructure",IsFpGroup);
DeclareSynonym("IsVHGroup",IsFpGroup and HasVHStructure);
DeclareAttribute("VerticalAction",IsVHGroup);
DeclareAttribute("HorizontalAction",IsVHGroup);
DeclareAttribute("VHRws",IsVHGroup);
DeclareGlobalFunction("VHGroup");
DeclareProperty("IsIrreducibleVHGroup",IsGroup);
DeclareAttribute("SimpleSubgroupVHGroup",IsVHGroup);
fi;

FR_LOCAL.VHSTRUCTURE := function(result,r,v,h)
    local i, m, n, getv, geth, addrel;
    
    m := Length(v);
    n := Length(h);
    getv := function(x)
        if x<0 then return 2*m+1-Position(v,-x); else
            return Position(v,x);
        fi;
    end;
    geth := function(x)
        if x<0 then return 2*n+1-Position(h,-x); else
            return Position(h,x);
        fi;
    end;
    result.trans := List([1..2*m],i->[]);
    result.out := List([1..2*m],i->[]);
    addrel := function(a,b,c,d)
        if IsBound(result.trans[a][b]) and (result.trans[a][b]<>c or result.out[a][b]<>d) then
            return true;
        fi;
        result.trans[a][b] := c;
        result.out[a][b] := d;
        return false;
    end;
    for i in r do
        if addrel(getv(i[1]),geth(-i[4]),getv(-i[3]),geth(i[2])) then
            return fail;
        fi;
        if addrel(getv(-i[1]),geth(i[2]),getv(i[3]),geth(-i[4])) then
            return fail;
        fi;
        if addrel(getv(i[3]),geth(-i[2]),getv(-i[1]),geth(i[4])) then
            return fail;
        fi;
        if addrel(getv(-i[3]),geth(i[4]),getv(i[1]),geth(-i[2])) then
            return fail;
        fi;
    od;
    if Set(result.trans,Length)<>[2*n] or Set(result.out,Length)<>[2*n] then
        return fail;
    fi;
    return true;
end;

InstallMethod(VHStructure, "for a f.p. group",
        [IsFpGroup],
        function(G)
    local i, v, h, r, result;
    
    v := [];
    h := [];
    r := [];
    for i in RelatorsOfFpGroup(G) do
        i := LetterRepAssocWord(i);
        if Length(i)<>4 then TryNextMethod(); fi;
        if AbsInt(i[2]) in v then
            i := i{[2,3,4,1]};
        fi;
        AddSet(v,AbsInt(i[1]));
        AddSet(h,AbsInt(i[2]));
        AddSet(v,AbsInt(i[3]));
        AddSet(h,AbsInt(i[4]));
        Add(r,i);
    od;
    if Intersection(v,h)<>[] then TryNextMethod(); fi;
    result := rec(v := GeneratorsOfGroup(G){v},
                  h := GeneratorsOfGroup(G){h});
    if FR_LOCAL.VHSTRUCTURE(result,r,v,h)=fail then
        TryNextMethod();
    fi;
    SetVHStructure(FamilyObj(One(G)),result);
    SetReducedMultiplication(G);
    return result;
end);

InstallMethod(VerticalAction, "for a VH group",
        [IsVHGroup],
        function(G)
    local r, m;
    r := VHStructure(G);
    m := MealyMachine(r.trans,r.out);
    SetAlphabetInvolution(m,[2*Length(r.h),2*Length(r.h)-1..1]);
    return GroupHomomorphismByImagesNC(Subgroup(G,r.v),SCGroup(m),
                   r.v,List([1..Length(r.v)],x->FRElement(m,x)));
end);

InstallMethod(HorizontalAction, "for a VH group",
        [IsVHGroup],
        function(G)
    local r, m;
    r := VHStructure(G);
    m := MealyMachine(TransposedMat(r.out),TransposedMat(r.trans));
    SetAlphabetInvolution(m,[2*Length(r.v),2*Length(r.v)-1..1]);
    return GroupHomomorphismByImagesNC(Subgroup(G,r.h),SCGroup(m),
                   r.h,List([1..Length(r.h)],x->FRElement(m,x)));
end);

InstallMethod(ViewObj, "for a VH group",
        [IsVHGroup], 10,
        function(G)
    local s, t;
    s := String(VHStructure(G).v);
    t := String(VHStructure(G).h);
    Print("<VH group on the generators ",s{[1..Length(s)-1]},"|",t{[2..Length(t)]},">");
end);

InstallMethod(FpElmKBRWS, "for a VH group",
        [IsElementOfFpGroupFamily and HasVHStructure],
        function(f)
    local r, iso, id, k;
    r := VHStructure(f);
    iso := IsomorphismFpMonoid(CollectionsFamily(f)!.wholeGroup);
    id := UnderlyingElement(Image(iso,One(f)));
    k := ShallowCopy(r.v);
    Append(k,List(r.v,x->ElementOfFpGroup(f,Inverse(UnderlyingElement(x)))));
    Append(k,r.h);
    Append(k,List(r.h,x->ElementOfFpGroup(f,Inverse(UnderlyingElement(x)))));
    k := ReducedConfluentRewritingSystem(Range(iso),ShortLexOrdering(FamilyObj(id),List(k,x->UnderlyingElement(Image(iso,x)))));
    return [iso,k,id];
end);

InstallMethod(FpElmEqualityMethod, "for a VH group",
        [IsElementOfFpGroupFamily and HasVHStructure],
        function(f)
    local iso,k,id;
    id:=FpElmKBRWS(f);
    iso:=id[1];k:=id[2];id:=id[3];
    return function(left,right)
        return ReducedForm(k,Gpword2MSword(id,UnderlyingElement(left),0))
               =ReducedForm(k,Gpword2MSword(id,UnderlyingElement(right),0));
    end;
end);

InstallMethod(FpElmComparisonMethod, "for a VH group",
        [IsElementOfFpGroupFamily and HasVHStructure],
        function(f)
    local iso,k,id;
    id:=FpElmKBRWS(f);
    iso:=id[1];k:=id[2];id:=id[3];
    return function(left,right)
        return ReducedForm(k,Gpword2MSword(id,UnderlyingElement(left),0))
               <ReducedForm(k,Gpword2MSword(id,UnderlyingElement(right),0));
    end;
end);
          
InstallOtherMethod(ReducedForm, "for a f.p. group element",
        [IsElementOfFpGroup],
        function(e)
    local f, k, x, y;
    f := FamilyObj(e);
    if IsBound(f!.reduce) and f!.reduce=true then
        k := FpElmKBRWS(f);
        x := Gpword2MSword(k[3],UnderlyingElement(e),0);
        y := ReducedForm(k[2],x);
        if x=y then return e; fi;
        return ElementOfFpGroup(f,MSword2gpword(UnderlyingElement(One(f)),y,0));
    else
        return e;
    fi;
end);

InstallMethod(InverseOp, "for an element of an f.p. group",
        [IsElementOfFpGroup],
        100,
        function(e)
    local f, k;
    f := FamilyObj(e);
    if IsBound(f!.reduce) and f!.reduce=true then
        k := FpElmKBRWS(f);
        return ElementOfFpGroup(f,MSword2gpword(UnderlyingElement(One(f)),ReducedForm(k[2],Gpword2MSword(k[3],Inverse(UnderlyingElement(e)),0)),0));
    else
        return ElementOfFpGroup(f,Inverse(UnderlyingElement(e)));
    fi;
end);

InstallGlobalFunction(VHGroup, function(arg)
    local l, i, m, n, v, h, r, f, addset;
    if Length(arg)=1 and IsList(arg[1]) then
        l := arg[1];
    else
        l := arg;
    fi;
    m := Maximum(List(l,x->Maximum(AbsInt(x[1]),AbsInt(x[3]))));
    n := Maximum(List(l,x->Maximum(AbsInt(x[2]),AbsInt(x[4]))));
    r := [];
    addset := function(p)
        if p in r then
            Error("Corner ",p," occurs too many times");
        fi;
        AddSet(r,p);
    end;
    for i in l do
        if Length(i)<>4 then
            Error("Bad length of relator ",i);
        fi;
        addset([i[1],i[2]]);
        addset([i[3],i[4]]);
        addset([-i[1],-i[4]]);
        addset([-i[3],-i[2]]);
    od;
    if Length(l)<>m*n or ForAny(l,x->0 in x) then
        Error("Missing corners ",Difference(Cartesian(Concatenation([-m..-1],[1..m]),Concatenation([-n..-1],[1..n])),r));
    fi;
    v := List([1..m],i->Concatenation("a",String(i)));
    h := List([1..n],i->Concatenation("b",String(i)));
    f := FreeGroup(Concatenation(v,h));
    v := GeneratorsOfGroup(f){[1..m]};
    h := GeneratorsOfGroup(f){[m+1..m+n]};
    f := f / List(l,x->v[AbsInt(x[1])]^SignInt(x[1])*h[AbsInt(x[2])]^SignInt(x[2])*v[AbsInt(x[3])]^SignInt(x[3])*h[AbsInt(x[4])]^SignInt(x[4]));
    i := rec(v := GeneratorsOfGroup(f){[1..m]},
             h := GeneratorsOfGroup(f){[m+1..m+n]});
    FR_LOCAL.VHSTRUCTURE(i,l,[1..m],[1..n]);
    SetVHStructure(f,i);
    SetVHStructure(FamilyObj(One(f)),i);
    SetReducedMultiplication(f);
    return f;
end);

RattaggiExample := rec(2_2 := VHGroup([1,1,-1,-1],[1,2,-1,-3],[1,3,2,-2],
                           [1,-3,-3,2],[2,1,-3,-2],[2,2,-3,-3],
                           [2,3,-3,1],[2,-3,3,2],[2,-1,-3,-1]),
                       # THM 2.3: (A6,A6), just infinite, irreducible
                       # CONJ 2.5: G0 is simple                       
                       2_15 := VHGroup([1,1,-1,-2],[1,2,-2,1],[1,3,-1,3],
                               [1,-2,2,-1],[2,1,-3,-3],[2,2,-3,3],
                               [2,3,-3,2],[2,-3,-3,1],[3,1,3,2]),
                       # THM 2.16: (A6,A6), just infinite, irreducible
                       # CONJ 2.17: G'' is simple, of index 192
                       2_18 := VHGroup([1,1,-2,-2],[1,2,-1,-1],[1,3,-1,-3],
                               [1,4,-1,-4],[1,5,-1,-6],[1,6,-1,-5],
                               [1,-1,2,2],[2,1,2,-3],[2,3,2,-4],
                               [2,4,-3,-5],[2,5,2,6],[2,-6,2,-2],
                               [2,-5,3,4],[3,1,-3,-2],[3,2,-3,-1],
                               [3,3,3,-6],[3,5,-3,-4],[3,6,3,-3]),
                       # THM 2.19: (A6,M12), just infinite, irreducible
                       # CONJ 2.20: G0 is simple                       
                       2_21 := VHGroup([1,1,-1,-1],[1,2,-1,-2],[1,3,-1,-4],
                               [1,4,-2,-3],[1,-4,-2,3],[2,1,-2,-2],
                               [2,2,-3,1],[2,3,-2,4],[2,-2,3,-1],
                               [3,1,3,-3],[3,2,3,-4],[3,3,3,4]),
                       # THM 2.22: (A6,S8), just infinite, irreducible
                       # CONJ 2.23: G0 is simple
                       2_26 := VHGroup([1,1,-1,-1],[1,2,-2,-3],[1,3,-1,-4],
                               [1,4,-1,-5],[1,5,-1,-6],[1,6,-1,-2],
                               [1,-2,2,3],[2,1,-2,-5],[2,2,2,-3],
                               [2,4,-2,4],[2,5,-2,-1],[2,6,-2,6]),
                       # THM 2.27: (A4,PSL(2,5)), irreducible, Lambda_2<>1, not residually finite
                       2_30 := VHGroup([1,1,-1,-1],[1,2,-2,-3],[1,3,-1,-4],
                               [1,4,-1,-5],[1,5,-1,-6],[1,6,-1,-2],
                               [1,7,2,-8],[1,8,2,8],[1,-8,2,-7],
                               [1,-7,3,7],[1,-2,2,3],[2,1,-2,-5],
                               [2,2,2,-3],[2,4,-2,4],[2,5,-2,-1],
                               [2,6,-2,6],[2,7,3,-7],[3,1,-3,8],
                               [3,2,-3,2],[3,3,-3,-4],[3,4,-3,1],
                               [3,5,-3,3],[3,6,-3,6],[3,8,-3,5]),
                       # THM 2.31: (A6,A16), virtually simple
                       2_33 := VHGroup([1,1,-1,-1],[1,2,-2,-3],[1,3,-1,-4],
                               [1,4,-1,-5],[1,5,-1,-6],[1,6,-1,-2],
                               [1,7,-2,-7],[1,-7,3,7],[1,-2,2,3],
                               [2,1,-2,-5],[2,2,2,-3],[2,4,-2,4],
                               [2,5,-2,-1],[2,6,-2,6],[2,7,-4,-7],
                               [3,1,4,4],[3,2,-3,-3],[3,3,-4,-2],
                               [3,4,4,7],[3,5,4,-6],[3,6,4,-1],
                               [3,-7,4,1],[3,-6,4,5],[3,-5,4,6],
                               [3,-4,4,-5],[3,-3,4,2],[3,-1,4,-4],
                               [4,3,4,-2]),
                       # THM 2.34: (ASL(3,2),A14), virtually simple
                       # CONJ 2.35: G0 is simple
                       2_36 := VHGroup([1,2,-1,-1],[2,2,-2,-1],[1,3,-2,-3],
                               [1,1,-2,-2],[2,1,-1,-3],[2,3,-1,-2]),
                       # THM 2.37: irreducible, not <b1,b2,b3>-separable
                       2_39 := VHGroup([1,2,-1,-1],[2,2,-2,-1],[1,3,-2,-3],
                               [1,1,-2,-2],[2,1,-1,-3],[2,3,-1,-2],
                               [3,2,-3,-1],[4,2,-4,-1],[3,3,-4,-3],
                               [3,1,-4,-2],[4,1,-3,-3],[4,3,-3,-2]),
                       # THM 2.40: a2/a1*a3/a4 \in N for all finite-index N
                       2_43 := VHGroup([1,1,-2,-2],[1,2,-1,-1],[1,3,-2,-3],
                               [1,4,2,-5],[1,5,-5,4],[1,-5,3,-4],
                               [1,-4,3,5],[1,-3,-2,2],[1,-1,-2,3],
                               [2,2,-2,-1],[2,4,-2,5],[2,5,4,-4],
                               [3,1,-4,-2],[3,2,-3,-1],[3,3,-4,-3],
                               [3,4,4,5],[3,-5,4,4],[3,-3,-4,2],
                               [3,-1,-4,3],[4,2,-4,-1],[4,-5,-5,-4],
                               [5,1,-5,3],[5,2,-5,-5],[5,3,-5,-1],
                               [5,4,-5,-2]),
                       # THM 2.44: (A10,A10), Z(a5)=<a5>, Z(a5^4) \ni b1
                       # THM 2.45: simple subgroup of index 4
                       2_46 := VHGroup([1,1,-2,-2],[1,2,-1,-1],[1,3,-2,-3],
                               [1,4,3,4],[1,-4,2,-4],[1,-3,-2,2],
                               [1,-1,-2,3],[2,2,-2,-1],[2,4,5,4],
                               [3,1,-4,-2],[3,2,-3,-1],[3,3,-4,-3],
                               [3,-4,-4,-4],[3,-3,-4,2],[3,-1,-4,3],
                               [4,2,-4,-1],[4,-4,5,-4],[5,1,-6,2],
                               [5,2,-6,-2],[5,3,-5,-3],[5,-2,-6,-1],
                               [5,-1,-6,1],[6,3,-6,-4],[6,4,-6,3]),
                       # THM 2.47: (M12,A8), G0 is simple
                       2_48 := VHGroup([1,1,-2,-2],[1,2,-1,-1],[1,3,-2,-3],
                               [1,4,2,-4],[1,5,2,-5],[1,6,-4,4],
                               [1,-6,4,6],[1,-5,-2,5],[1,-4,-4,-6],
                               [1,-3,-2,2],[1,-1,-2,3],[2,2,-2,-1],
                               [2,4,-3,-6],[2,6,-3,-4],[2,-6,3,6],
                               [3,1,-4,-2],[3,2,-3,-1],[3,3,-4,-3],
                               [3,4,5,5],[3,5,-4,-4],[3,-5,-4,-5],
                               [3,-3,-4,2],[3,-1,-4,3],[4,2,-4,-1],
                               [4,-4,5,-5],[5,1,-5,-1],[5,2,-5,2],
                               [5,3,-5,5],[5,4,-5,-3],[5,6,-5,6]),
                       # THM 2.49: (A10,A12), simple subgroup of index 12
                       2_50 := VHGroup([1,1,-2,-2],[1,2,-1,-1],[1,3,-2,-3],
                               [1,4,3,4],[1,5,-1,-5],[1,-4,2,-4],
                               [1,-3,-2,2],[1,-1,-2,3],[2,2,-2,-1],
                               [2,4,4,4],[2,5,-5,-5],[2,-5,-5,5],
                               [3,1,-4,-2],[3,2,-3,-1],[3,3,-4,-3],
                               [3,5,4,-4],[3,-5,4,-5],[3,-4,4,5],
                               [3,-3,-4,2],[3,-1,-4,3],[4,2,-4,-1],
                               [5,1,-5,-3],[5,2,-5,-2],[5,3,-5,4],
                               [5,4,-5,1]),
                       # THM 2.51: (A10,10), simple subgroup of index 40
                       2_52 := VHGroup([1,1,-2,-2],[1,2,-1,-1],[1,3,-2,-3],
                               [1,4,1,5],[1,-5,2,-5],[1,-4,-4,-4],
                               [1,-3,-2,2],[1,-1,-2,3],[2,2,-2,-1],
                               [2,4,2,5],[2,-4,-3,-4],[3,1,-4,-2],
                               [3,2,-3,-1],[3,3,-4,-3],[3,5,4,-4],
                               [3,-5,-5,-5],[3,-4,4,5],[3,-3,-4,2],
                               [3,-1,-4,3],[4,2,-4,-1],[4,-5,5,-5],
                               [5,1,5,4],[5,2,-5,3],[5,3,-5,2],
                               [5,-4,5,-1]),
                       # PROP 2.53: (3840,S10), not residually finite, irreducible
                       # THM 2.54: G0 has no f.i. subgroup, not simple
                       2_56 := VHGroup([1,1,-2,-2],[1,2,-1,-1],[1,3,-2,-3],
                               [1,4,-2,4],[1,-4,-2,-4],[1,-3,-2,2],
                               [1,-1,-2,3],[2,2,-2,-1],[3,1,-4,-2],
                               [3,2,-3,-1],[3,3,-4,-3],[3,4,-3,4],
                               [3,-3,-4,2],[3,-1,-4,3],[4,2,-4,-1],
                               [4,4,-4,-4]),
                       # THM 2.57: if w=a2/a1*a3/a4, then G/<w^2> is not residually finite, and not virtually torsion-free.
                       2_58 := VHGroup([1,1,-1,2],[1,2,-2,-3],[1,3,-2,1],
                               [1,4,-2,-5],[1,5,-2,5],[1,-5,-2,-4],
                               [1,-4,2,-1],[1,-3,-2,3],[1,-2,2,4],
                               [2,1,-3,2],[2,2,-3,1],[3,1,3,2],
                               [3,3,-3,-3],[3,4,3,-4],[3,5,-3,5]),
                       # THM 2.59: (A6,S5), SQ-universal, irreducible
                       # CONJ 2.61: intersection of all N = G0
                       # CONJ 2.63: QZ(H2)=1
                       # CONJ 2.65: G/N has (T), for infinite-index non-trivial N
                       2_70 := VHGroup([1,1,-1,-2],[1,2,-2,-1],[1,3,-2,1],
                               [1,-3,2,3],[1,-2,-2,-3],[2,1,-2,2])
                       # CONJ 2.70: G0 is simple
                       );
