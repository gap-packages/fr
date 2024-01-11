#############################################################################
##
#W helpers.gi                                               Laurent Bartholdi
##
#Y Copyright (C) 2012-2016, Laurent Bartholdi
##
#############################################################################
##
##  This file contains helper code for functionally recursive groups,
##  in particular related to the geometry of groups.
##
#############################################################################

BindGlobal("INSTALLPRINTERS@", function(filter)
    InstallMethod(PrintObj, "(FR)", [filter], 2*SUM_FLAGS, function(x) Print(String(x)); end);
    InstallMethod(ViewObj, "(FR)", [filter], 2*SUM_FLAGS, function(x) Print(ViewString(x)); end);
    InstallMethod(Display, "(FR)", [filter], 2*SUM_FLAGS, function(x) Print(DisplayString(x)); end);
end);
#############################################################################

#############################################################################
##
#F Products
##
InstallGlobalFunction(TensorSum, function(arg)
    local d;
    if Length(arg) = 0 then
        Error("<arg> must be nonempty");
    elif Length(arg) = 1 and IsList(arg[1])  then
        if arg[1]=[]  then
            Error("<arg>[1] must be nonempty");
        fi;
        arg := arg[1];
    fi;
    d := TensorSumOp(arg,arg[1]);
    if ForAll(arg, HasSize) then
        if ForAll(arg, IsFinite) then
            SetSize(d, Product( List(arg, Size)));
        else
            SetSize(d, infinity);
        fi;
    fi;
    return d;
end);
#############################################################################

#############################################################################
##
#H WordGrowth(g,options)
##
COLOURLIST@ :=
  ["red","blue","green","gray","yellow","cyan","orange","purple"];
BindGlobal("COLOURS@", function(i)
    return COLOURLIST@[(i-1) mod Length(COLOURLIST@)+1];
end);

BindGlobal("EXEC@", rec());
# If 1 argument is passed, search path to find appropriate executable.
# If >1 arguments, the first is a generic name, such as "psviewer", and the
#     other arguments are possibilities to be searched in sequence, in the
#     form ["progname","arg1",...,"argn"]. If the last argument is true, then
#     do not raise an error if nothing was found.
BindGlobal("CHECKEXEC@", function(arg)
    local a, prog, command;
    prog := arg[1];
    
    if IsBound(EXEC@.(prog)) then return; fi;
    
    if Length(arg)=1 then
        command := IO_FindExecutable(prog);
    else
        for a in arg{[2..Length(arg)]} do
            if a=true then return; fi;
            command := IO_FindExecutable(a[1]);
            if command<>fail then
                command := JoinStringsWithSeparator(Concatenation([command],a{[2..Length(a)]})," ");
                break;
            fi;
        od;
    fi;
    
    while command=fail do
        Error("Could not find program \"",prog,"\" -- make sure it is installed, and/or set manually EXEC@FR.",prog);
    od;
    EXEC@.(prog) := command;
end);

BindGlobal("OUTPUTTEXTSTRING@", function(s)
    local f;
    f := OutputTextString(s,false);
    SetPrintFormattingStatus(f,false);
    return f;
end);

BindGlobal("STRINGGROUP@",
        function(O)
    local s, os;
    s := "";
    os := OutputTextString(s,true);
    PrintTo(os,O);
    CloseStream(os);
    return s;
end);

BindGlobal("EXECINSHELL@", function(input,command,detach)
    local tmp, outs, output;
    outs := "";
    output := OUTPUTTEXTSTRING@(outs);
    CHECKEXEC@("sh");

    if detach=fail then
        if IsString(input) then input := InputTextString(input); fi;
    else
	tmp := Filename(DirectoryTemporary(), "stdin");
 
        if not IsString(input) then
	    input := ReadAll(input);
        fi;
        WriteAll(OutputTextFile(tmp,false), input);
        input := InputTextNone();
        CHECKEXEC@("cat");
        command := Concatenation("cat ",tmp,"|",command,"&");
    fi;
    Process(DirectoryCurrent(), EXEC@.sh, input, output, ["-c", command]);
    return outs;
end);

BindGlobal("DOT2DISPLAY@", function(str,prog)
    local command;
    
    CHECKEXEC@(prog);
    CHECKEXEC@("psviewer",["display","-flatten","-"],["gv","-"],true);
    if ValueOption("usesvg")<>fail or EXEC@.psviewer=fail then
        CHECKEXEC@("svgviewer",["svg-view","--stdin"]);
        command := Concatenation(EXEC@.(prog)," -Tsvg 2>/dev/null | ",EXEC@.svgviewer);
    else
        command := Concatenation(EXEC@.(prog)," -Gbgcolor=white -Tps 2>/dev/null | ",EXEC@.psviewer);
    fi;
    return EXECINSHELL@(str,command,ValueOption("detach"));
end);

BindGlobal("APPEND@", function(arg)
    local i;
    for i in [2..Length(arg)] do
        Append(arg[1],String(arg[i]));
    od;
end);

BindGlobal("CONCAT@", function(arg)
    local i, s;
    s := "";
    for i in arg do
        Append(s,String(i));
    od;
    return s;
end);

InstallGlobalFunction(WordGrowth, function(arg)
    local gpgens, gens, sphere, i, j, k, n, t, result, limit, keep, g, act,
          tile, point, draw, S, plotedge, plotvertex,
          trackgroup, trackgens, track, trackhom,
          group, options, optionnames;
    
    optionnames := ["track","limit","draw","point","ball","sphere","balls",
                    "spheres","spheresizes","ballsizes","act","tile"];
    if Length(arg)=2 then
        group := arg[1];
        options := arg[2];
    elif Length(arg)>2 then
        Error("Too many arguments for WordGrowth");
    else
        group := arg[1];
        options := rec();
        for i in optionnames do
            if ValueOption(i)<>fail then
                options.(i) := ValueOption(i);
            fi;
        od;
    fi;
    
    if IsInt(options) then
        options := rec(spheresizes := options);
    fi;
    i := Difference(RecNames(options),optionnames);
    if i<>[] then
        Info(InfoFR,1,"WordGrowth: unused options ",i);
    fi;

    if IsGroup(group) then
        g := group;
        gpgens := Set(GeneratorsOfGroup(g));
        gens := Union(gpgens,List(gpgens,Inverse),[One(g)]);
    elif IsSemigroup(group) then
        g := group;
        gens := Set(GeneratorsOfSemigroup(g));
    elif IsList(group) then
        g := Semigroup(group);
        gens := Set(group);
    else
        TryNextMethod();
    fi;

    if IsBound(options.track) then
        track := [];
        if IsGroup(g) and not IsList(group) then
            if IsList(options.track) then
                trackgroup := FreeGroup(options.track);
            else
                trackgroup := FreeGroup(Length(GeneratorsOfGroup(g)));
            fi;
            trackgens := GroupHomomorphismByImages(trackgroup,g,GeneratorsOfGroup(trackgroup),GeneratorsOfGroup(g));
            trackgens := List(gens,x->PreImagesRepresentativeNC(trackgens,x));
            trackhom := GroupHomomorphismByImagesNC(trackgroup,g,GeneratorsOfGroup(trackgroup),GeneratorsOfGroup(g));
        elif IsMonoid(g) and not IsList(group) then
            if IsList(options.track) then
                trackgroup := FreeMonoid(options.track);
            else
                trackgroup := FreeMonoid(Length(GeneratorsOfMonoid(g)));
            fi;
            trackgens := GeneratorsOfMonoid(trackgroup){List(gens,x->Position(GeneratorsOfMonoid(g),x))};
            trackhom := FreeMonoidNatHomByGeneratorsNC(trackgroup,g);
        elif IsSemigroup(g) then
            if IsList(options.track) then
                trackgroup := FreeSemigroup(options.track);
            else
                trackgroup := FreeSemigroup(Length(GeneratorsOfSemigroup(g)));
            fi;
            trackgens := GeneratorsOfSemigroup(trackgroup){List(gens,x->Position(GeneratorsOfSemigroup(g),x))};
            trackhom := FreeSemigroupNatHomByGeneratorsNC(trackgroup,g);
        fi;
    else
        track := fail;
    fi;

    if IsBound(options.limit) then
        limit := options.limit;
    else
        limit := infinity;
    fi;
    keep := IsBound(options.ball) or IsBound(options.balls) or
            IsBound(options.spheres) or not IsBound(gpgens);
    if IsBound(options.point) then
        point := options.point;
    else
        point := fail;
    fi;
    if IsBound(options.draw) then
        draw := options.draw;
    else
        draw := fail;
    fi;
    if IsBound(options.sphere) and limit=infinity then
        limit := options.sphere;
    fi;
    if IsBound(options.spheres) and limit=infinity then
        limit := options.spheres;
    fi;
    if IsBound(options.spheresizes) and limit=infinity then
        limit := options.spheresizes;
    fi;
    if IsBound(options.ball) and limit=infinity then
        limit := options.ball;
    fi;
    if IsBound(options.balls) and limit=infinity then
        limit := options.balls;
    fi;
    if IsBound(options.ballsizes) and limit=infinity then
        limit := options.ballsizes;
    fi;
    if IsBound(options.act) then
        act := options.act;
    elif point=fail then
        act := OnRight;
    else
        act := OnPoints;
    fi;
    if IsBound(options.tile) then
        tile := options.tile;
    else
        tile := false;
    fi;
    if draw<>fail then
        S := "digraph cayley {\n";
        plotedge := function(nsrc,src,ndst,dst,gen)
            local dir, col;
            dir := "forward";
            if IsBound(gpgens) then
                if IsOne(gen^2) then
                    dir := "both";
                    if ndst=nsrc and dst<src then return; fi;
                fi;
                if gen in gpgens then
                    col := Position(gpgens,gen);
                else
                    col := Position(gpgens,gen^-1);
                    dir := "back";
                fi;
            else
                col := Position(gens,gen);
            fi;
            APPEND@(S,"  ",nsrc,".",src," -> ",ndst,".",dst," [color=",COLOURS@(col),",dir=",dir,"];\n");
        end;
        plotvertex := function(nsrc,src)
            APPEND@(S,"  ",nsrc,".",src," [height=0.3,width=0.6,fixedsize=true]\n");
        end;
    fi;

    i := PositionProperty(gens,IsMultiplicativeElementWithOne and IsOne);
    if i<>fail then
        sphere := [gens[i]];
        if track<>fail then
            Add(track,[trackgens[i]]);
        fi;
        Remove(gens,i);
    elif HasOne(g) then
        sphere := [One(g)];
        if track<>fail then
            Add(track,[One(trackgroup)]);
        fi;
    else
        sphere := [];
        if track<>fail then Add(track,[]); fi;
    fi;
    if point=fail then
        sphere := [sphere,Difference(gens,sphere)];
        if track<>fail then
            Add(track,trackgens);
            if track[1]<>[] then
                t := Position(track[2],track[1][1]);
                if t<>fail then Remove(track[2],t); fi;
            fi;
        fi;
    else
        sphere := DifferenceLists(List(sphere,g->act(point,g)),[fail]);
        if track=fail then
            sphere := [sphere,[]];
            for i in gens do
                j := act(point,i);
                if j=fail then continue; fi;
                if not j in sphere[1] then AddSet(sphere[2],j); fi;
            od;
        else
            sphere := [sphere,[]];
            Add(track,[]);
            for i in [1..Length(gens)] do
                j := act(point,gens[i]);
                if j=fail then continue; fi;
                if not (j in sphere[1] or j in sphere[2]) then
                    t := PositionSorted(sphere[2],j);
                    Add(sphere[2],j,t);
                    Add(track,trackgens[i],t);
                fi;
            od;
        fi;
    fi;
    if draw<>fail then
        if sphere[1]<>[] then
            for n in [1..Length(gens)] do
                if point=fail then
                    plotedge(0,1,1,n,gens[n]);
                else
                    j := act(point,gens[n]);
                    i := Position(sphere[2],j);
                    if i<>fail then
                        plotedge(0,1,1,i,gens[n]);
                    elif j=point then
                        plotedge(0,1,0,1,gens[n]);
                    fi;
                fi;
            od;
        fi;
        if limit<infinity then limit := limit+1; fi;
    fi;
    result := List(sphere,Size);

    n := 1; while n < limit do
        Add(sphere,[]);
        if track<>fail then
            Add(track,[]);
        fi;
        if track<>fail then
            for i in [1..Length(gens)] do for j in [1..Length(sphere[n+1])] do
                k := act(sphere[n+1][j],gens[i]);
                if k=fail then continue; fi;
                if (IsBound(gpgens) and
                    not (k in sphere[n] or k in sphere[n+1])) or
                   (not IsBound(gpgens) and not ForAny([1..n+1],i->k in sphere[i])) then
                    t := PositionSorted(sphere[n+2],k);
                    if not IsBound(sphere[n+2][t]) or sphere[n+2][t]<>k then
                        Add(sphere[n+2],k,t);
                        Add(track[n+2],track[n+1][j]*trackgens[i],t);
                    fi;
                fi;
            od; od;
        elif draw=fail then
            for i in gens do for j in sphere[n+1] do
                k := act(j,i);
                if (IsBound(gpgens) and
                    not (k in sphere[n] or k in sphere[n+1])) or
                   (not IsBound(gpgens) and not ForAny([1..n+1],i->k in sphere[i])) then
                    AddSet(sphere[n+2],k);
                fi;
            od; od;
        else
            for i in gens do for j in [1..Length(sphere[n+1])] do
                k := act(sphere[n+1][j],i);
                if k=fail then continue; fi;
                t := 1; while t <= Length(sphere) do
                    if IsBound(sphere[t]) and k in sphere[t] then break; fi;
                    t := t+1;
                od;
                if t>Length(sphere) then
                    if n+1<limit then
                        Add(sphere[n+2],k); t := n+2;
                    else
                        continue;
                    fi;
                fi;
                if (not IsBound(gpgens)) or ((i in gpgens and t >= n+1) or (t >= n+2)) then
                    plotedge(n,j,t-1,Position(sphere[t],k),i);
                fi;
            od; od;
        fi;
        if limit=infinity and sphere[n+2]=[] then
            Remove(sphere);
            if sphere[n+1]=[] then Remove(sphere); Remove(result); fi;
            if track<>fail then
                Remove(track);
                if track[n+1]=[] then Remove(track); fi;
            fi;
            break;
        fi;
        Add(result,Size(sphere[n+2]));
        n := n+1;
        if not keep then Unbind(sphere[n-1]); fi;
    od;
    if limit=0 then
        Unbind(result[2]); Unbind(sphere[2]);
    fi;
    if IsBound(options.spheresizes) then
        return result;
    elif IsBound(options.ballsizes) then
        return List([1..Length(result)],i->Sum(result{[1..i]}));
    elif IsBound(options.sphere) then
        if track=fail then
            return sphere[limit+1];
        else
            return [sphere[limit+1],trackhom,track[limit+1]];
        fi;
    elif IsBound(options.spheres) then
        if track=fail then return sphere; else return [sphere,trackhom,track]; fi;
    elif IsBound(options.ball) then
        if track=fail then
            return Union(sphere);
        else
            sphere := Concatenation(sphere); track := Concatenation(track);
            SortParallel(sphere,track);
            return [sphere,trackhom,track];
        fi;
    elif IsBound(options.balls) then
        if track=fail then
            return List([1..Length(sphere)],i->Union(sphere{[1..i]}));
        else
            t := [[],trackhom,[]];
            for i in [1..Length(sphere)] do
                Add(t[1],Concatenation(sphere{[1..i]}));
                Add(t[3],Concatenation(track{[1..i]}));
                SortParallel(t[1][i],t[3][i]);
            od;
            return t;
        fi;
    elif IsBound(options.indet) then
        i := options.indet;
        if HasIsUnivariatePolynomial(i) and IsUnivariatePolynomial(i) then
            i := IndeterminateNumberOfUnivariateRationalFunction(i);
        else i := 1; fi;
        return UnivariatePolynomial(Integers,result,i);
    elif draw<>fail then
        if not IsBound(result[n+1]) then n := n-1; fi;
        for n in [0..n] do for i in [1..result[n+1]] do
            plotvertex(n,i);
        od; od;
        Append(S,"}\n");
        if IsString(draw) then
            AppendTo(draw,S);
        else
            if IsBoundGlobal("JupyterRenderable") then
                return ValueGlobal("JupyterRenderable")(rec(("image/svg+xml") :=IO_PipeThrough("dot",["-Tsvg"],S)),rec());
            else
                DOT2DISPLAY@(S, "neato");
            fi;
        fi;
    else
        return result; # by default, same as 'spheresizes'
    fi;
end);

InstallOtherMethod(Draw, "(FR) default",
        [IsObject],
        function(l)
    if IsBoundGlobal("JupyterRenderable") then
        return WordGrowth(l,rec(draw:=true));
    else
        WordGrowth(l,rec(draw:=true));
    fi;
end);

InstallOtherMethod(Draw, "(FR) default, with filename",
        [IsObject,IsString],
        function(l,w)
    WordGrowth(l,rec(draw:=w));
end);

InstallOtherMethod(Draw, "(FR) default, with options",
        [IsObject,IsRecord],
        function(l,options)
    options := ShallowCopy(options);
    options.draw := true;
    if IsBoundGlobal("JupyterRenderable") then # a hack
        return WordGrowth(l,options);
    else
        WordGrowth(l,options);
    fi;
end);

InstallMethod(Ball, "(FR) for an object and a limit radius",
        [IsObject,IsInt],
        function(x,n)
    return WordGrowth(x,rec(ball:=n));
end);

InstallMethod(Sphere, "(FR) for an object and a limit radius",
        [IsObject,IsInt],
        function(x,n)
    return WordGrowth(x,rec(sphere:=n));
end);

InstallGlobalFunction(OrbitGrowth, # "(FR) for an object and point or options",
        function(arg)
    if Length(arg)=2 then
        return WordGrowth(arg[1],rec(point:=arg[2]));
    else
        return WordGrowth(arg[1],rec(point:=arg[2],limit:=arg[3]));
    fi;
end);
#############################################################################

#############################################################################
##
#H Draw order relations
##
BindGlobal("ORDER2DOT@", function(R)
    local i, j, succ, S;
    
    if not IsBinaryRelationOnPointsRep(R) then
        R := AsBinaryRelationOnPoints(R);
    fi;

    S := "digraph ";
    if HasName(R) and ForAll(Name(R),IsAlphaChar) then
        APPEND@(S, "\"",Name(R),"\"");
    else
        Append(S,"HasseDiagram");
    fi;
    Append(S," {\n");
    for i in [1..DegreeOfBinaryRelation(R)] do
        APPEND@(S,i," [shape=circle]\n");
    od;
    
    succ := Successors(R);

    for i in [1..DegreeOfBinaryRelation(R)] do
        for j in succ[i] do
            APPEND@(S,"  ",i," -> ",j," [label=\".\"];\n");
        od;
    od;
    Append(S,"}\n");
    return S;
end);

InstallMethod(Draw, "(FR) for a binary relation",
        [IsBinaryRelation],
        function(R)
    DOT2DISPLAY@(ORDER2DOT@(R),"dot");
end);

InstallMethod(Draw, "(FR) for a binary relation and a filename",
        [IsBinaryRelation,IsString],
        function(R,S)
    AppendTo(S,ORDER2DOT@(R));
end);

InstallMethod(HeightOfPoset, "(FR) for a binary relation",
        [IsBinaryRelation],
        function(poset)
  local s, n, min;
  s := Elements(Source(poset));
  n := -1;
  repeat
    min := Filtered(s,x->Intersection(ImagesElm(poset,x),s)=[x]);
    s := Difference(s,min);
    n := n+1;
  until s=[];
  return n;
end);
#############################################################################

#############################################################################
## forward orbits
InstallMethod(ForwardOrbit, "(FR) forward orbit under a group element",
	[IsMultiplicativeElementWithInverse,IsObject],
        function(g,x)
    local l, y;
    if Inverse(g)=fail then TryNextMethod(); fi;
    l := [];
    y := x;
    repeat
        Add(l,y);
        y := y^g;
    until y=x;
    return l;
end);

InstallMethod(ForwardOrbit, "(FR) forward orbit under a semigroup element",
	[IsMultiplicativeElement,IsObject],
        function(g,x)
    local l, n;
    l := [];
    n := 0;
    repeat
        Add(l,x);
        x := x^g;
        Add(l,x);
        x := x^g;
        n := n+1;
    until x=l[n];
    # now we have l of length 2n; and l[2n+1]=l[n]. make it minimal.
    l := PeriodicList(l{[1..n-1]},l{[n..2*n]});
    CompressPeriodicList(l);
    return Concatenation(PrePeriod(l),Period(l));
  end
);
#############################################################################

#############################################################################
##
#H StringByInt
#H Rename subobjects if they have "=" (but not IsIdentical) named objects
##
InstallGlobalFunction(StringByInt, function(arg)
    local base, result, digit, n;
    if Size(arg)=1 then base := 2; else base := arg[2]; fi;
    result := "";
    n := arg[1];
    while n>0 do
        digit := n mod base;
        if digit <= 10 then
            digit := CHAR_INT(digit+INT_CHAR('0'));
        else
            digit := CHAR_INT(digit+INT_CHAR('a')-10);
        fi;
        Add(result,digit,1);
        n := QuoInt(n,base);
    od;
    return result;
end);

InstallGlobalFunction(PositionInTower, function(seq,x)
    local low, high, mid;
    low := 1; high := Size(seq);
    if not IsList(x) then x := [x]; fi;
    if x=seq[Length(seq)] then
        return infinity;
    elif not IsSubset(seq[1],x) then
        return fail;
    fi;
    while low < high-1 do
        mid := QuoInt(low+high,2);
        if IsSubset(seq[mid],x) then
            low := mid;
        else high := mid; fi;
    od;
    return low;
end);

InstallMethod(RenameSubobjects, "(FR) for an object and a list of named objs",
        [IsObject, IsList],
        function(obj,refobj)
    local i;
    if IsList(obj) then
        for i in obj do RenameSubobjects(i,refobj); od;
    elif IsRecord(obj) then
        for i in RecNames(obj) do RenameSubobjects(obj.(i),refobj); od;
    elif not HasName(obj) then
        i := Position(refobj,obj);
        if i<>fail then SetName(obj,Name(refobj[i])); fi;
    fi;
end);

InstallGlobalFunction(CoefficientsInAbelianExtension, function(x,seq,G)
    local ord, i, j, k;
    ord := [];
    for i in [1..Length(seq)] do
        j := 1; k := seq[i];
        while not k in G do j := j+1; k := k*seq[i]; od;
        Add(ord,[0..j-1]);
    od;
    return Filtered(Cartesian(ord),
                   s->Product([1..Length(seq)],n->seq[n]^s[n])/x in G);
end);

BindGlobal("MAPPEDWORD@", function(arg)
    local i, e, w, gens;

    w := arg[1];
    while not IsAssocWord(w) do w := UnderlyingElement(w); od;
    w := LetterRepAssocWord(w);
    gens := arg[2];
    if w=[] then
        if Length(arg)=2 then return One(gens[1]); else return arg[3]; fi;
    fi;
    e := fail;
    for i in w do
        if e=fail then
            if i>0 then e := gens[i]; else e := Inverse(gens[-i]); fi;
        elif i>0 then
            e := e*gens[i];
        elif i<0 then
            e := e/gens[-i];
        fi;
    od;
    return e;
    #!!! could be a bit smarter here: if all gens are free group elements,
    # can construct a word by concatenation and convert it once to a free
    # group element; if all gens are FR elements on the same machine, idem.
    # this would presumably speed up quite a lot the code.
end);

InstallGlobalFunction(MagmaEndomorphismByImagesNC, function(f,im)
    return MagmaHomomorphismByImagesNC(f,f,im);
end);

InstallGlobalFunction(MagmaHomomorphismByImagesNC, function(f,g,im)
    local one;

    if IsGroup(f) and IsGroup(g) then
        return GroupHomomorphismByImagesNC(f,g,GeneratorsOfGroup(f),im);
    elif IsFreeGroup(f) or (HasIsFreeMonoid(f) and IsFreeMonoid(f)) or (HasIsFreeSemigroup(f) and IsFreeSemigroup(f)) then
        if IsMonoid(g) then
            one := One(g);
        else
            one := fail;
        fi;
        return MagmaHomomorphismByFunctionNC(f,g,w->MAPPEDWORD@(w,im,one));
    fi;

    if IsMonoid(f) and IsMonoid(g) then
        return MagmaHomomorphismByFunctionNC(f,g,function(x)
            local s;
            s := ShortMonoidWordInSet(f,x,infinity);
            if Length(s)<2 then return fail; fi;
            return MAPPEDWORD@(s[2],im,One(g));
        end);
    else
        return MagmaHomomorphismByFunctionNC(f,g,function(x)
            local s;
            s := ShortSemigroupWordInSet(f,x,infinity);
            if Length(s)<2 then return fail; fi;
            return MAPPEDWORD@(s[2],im,fail);
        end);
    fi;
end);

InstallMethod(ImagesSet, "(FR) for a magma homomorphism",
        [IsMagmaHomomorphism, IsMagma],
	function(map,set)
    local f;
    if HasGeneratorsOfGroup(set) then
        f := GeneratorsOfGroup;
    elif HasGeneratorsOfMonoid(set) or HasGeneratorsOfSemigroup(set) or HasGeneratorsOfMagma(set) then
        f := GeneratorsOfMagma;
    else
    	TryNextMethod();
    fi;
    if not IsGroupHomomorphism(map) then
        f := GeneratorsOfMagma;
    fi;
    set := List(f(set),x->x^map);
    if f=GeneratorsOfGroup then
        return GroupByGenerators(set);
    else
    	return MagmaByGenerators(set);
    fi;
end);
#############################################################################

#############################################################################
##
#H ShortMonoidRelations
##
BindGlobal("FINDMONOIDRELATIONS@", function(gens,n)
    local free, seen, reducible, i, iterate, result, freegens;
    free := FreeMonoid(Length(gens));
    freegens := GeneratorsOfMonoid(free);
    seen := NewDictionary(gens[1],true);
    reducible := NewDictionary(freegens[1],false);
    iterate := function(level,elem,word)
        local i, newword;
        if level=0 then
            if KnowsDictionary(seen,elem) then
                Add(result,[word,LookupDictionary(seen,elem)]);
                AddDictionary(reducible,word);
            fi;
            AddDictionary(seen,elem,word);
        else
            AddDictionary(seen,elem,word);
            for i in [1..Length(gens)] do
                newword := word*freegens[i];
                if ForAll([1..Length(newword)],i->not KnowsDictionary(reducible,Subword(newword,i,Length(newword)))) then
                    iterate(level-1,elem*gens[i],newword);
                fi;
            od;
        fi;
    end;
    result := [FreeMonoidNatHomByGeneratorsNC(free,MonoidByGenerators(gens))];
    for i in [1..n] do iterate(i,gens[1]^0,freegens[1]^0); od;
    return result;
end);

InstallMethod(ShortMonoidRelations, "for a monoid and a length",
        [IsMonoid,IsInt],
        function(m,n)
    return FINDMONOIDRELATIONS@(GeneratorsOfMonoid(m),n);
end);
InstallMethod(ShortMonoidRelations, "for a list and a length",
        [IsListOrCollection,IsInt],
        FINDMONOIDRELATIONS@);

BindGlobal("FINDGROUPRELATIONS_ADD@", function(result,w)
    if Sum(ExponentSums(w))>=0 then
        Add(result,w);
    else
        Add(result,w^-1);
    fi;
end);

BindGlobal("FINDGROUPRELATIONS@", function(group,gens,n)
    local freegens, pigens, freegensinv, pi, inv,
          i, j, k, m, mm, newf, newg, result, rws, seen, todo, x, y, z;

    gens := Unique(Concatenation(gens,List(gens,Inverse)));
    result := [EpimorphismFromFreeGroup(group)];
    y := GeneratorsOfGroup(group);
    z := GeneratorsOfGroup(Source(result[1]));
    for i in [1..Length(y)] do
        for j in [i+1..Length(y)] do
            if y[i]=y[j] then
                FINDGROUPRELATIONS_ADD@(result,z[i]/z[j]);
            elif y[i]=y[j]^-1 then
                FINDGROUPRELATIONS_ADD@(result,z[i]*z[j]);
            fi;
        od;
    od;
    x := FreeMonoid(Length(gens));
    rws := KnuthBendixRewritingSystem(x/[]);
    freegens := GeneratorsOfMonoid(x);
    pigens := List(gens,x->First(GeneratorsOfSemigroup(Source(result[1])),y->y^result[1]=x));
    pi := MagmaHomomorphismByImagesNC(x,Source(result[1]),pigens);
    freegensinv := [];
    for i in [1..Length(gens)] do
        j := Position(gens,gens[i]^-1);
        Add(freegensinv,freegens[j]);
        AddRuleReduced(rws,[[i,j],[]]);
        if j=i then
            FINDGROUPRELATIONS_ADD@(result,(freegens[i]^pi)^2);
        fi;
    od;
    freegensinv := List([1..Length(freegens)],
                        i->freegens[Position(gens,gens[i]^-1)]);
    inv := MagmaHomomorphismByFunctionNC(x,x,
                   w->Reversed(MappedWord(w,freegens,freegensinv)));
    seen := NewDictionary(gens[1],true);

    todo := NewFIFO([[0,freegens[1]^0,gens[1]^0]]); # store [len, normalform, elt]
    AddDictionary(seen,gens[1]^0,freegens[1]^0);

    for i in todo do
        m := i[1]+1;
        for j in [1..Length(gens)] do
            newf := i[2]*freegens[j];
            if not IsReducedForm(rws,newf) then continue; fi;
            newg := i[3]*gens[j];
            x := LookupDictionary(seen,newg);
            if x<>fail then
                mm := Length(x);
                newf := ReducedForm(rws,newf*x^inv);
                if Length(newf)<m+mm then continue; fi;
                FINDGROUPRELATIONS_ADD@(result,newf^pi);
                x := [newf^2,(newf^inv)^2];
                for j in [1..2] do
                    if m=mm then
                        for k in [1..2*m] do
                            y := ReducedForm(rws,Subword(x[j],k,k+m-1));
                            z := ReducedForm(rws,Subword(x[3-j],2*m-k+2,3*m-k+1));
                            if y>z then y := [y,z]; else y := [z,y]; fi;
                            AddRuleReduced(rws,List(y,LetterRepAssocWord));
                        od;
                    fi;
                    for k in [1..m+mm] do
                        y := Subword(x[j],k,k+mm);
                        z := Subword(x[3-j],m+mm-k+2,2*m+mm-k);
                        AddRuleReduced(rws,List([y,z],LetterRepAssocWord));
                    od;
                od;
            elif m<=n then
                AddDictionary(seen,newg,newf);
                Add(todo,[m,newf,newg]);
            fi;
        od;
    od;
    return result;
end);

InstallMethod(ShortGroupRelations, "for a group and a length",
        [IsGroup,IsInt],
        function(g,n)
    return FINDGROUPRELATIONS@(g,GeneratorsOfGroup(g),n);
end);
InstallMethod(ShortGroupRelations, "for a list and a length",
        [IsListOrCollection,IsInt],
        function(g,n)
    return FINDGROUPRELATIONS@(Group(g),g,n);
end);

BindGlobal("SHORTWORDINSET@", function(f,fgen,fone,ggen,gone,set,result,n)
    local x, newfx, newlen, i, seen, forbidden, todo, justone, compare;
    if n<0 then
        n := -n;
        justone := false;
    else
        justone := true;
    fi;
    compare := function(set,elm)
        if IsFunction(set) then
            return set(elm);
        elif FamilyObj(elm)=FamilyObj(set) and elm=set then
            return true;    
        elif IsListOrCollection(set) then
            return elm in set;
        fi;
        return false;
    end;
    if fone=fail then
        todo := NewFIFO(List([1..Length(fgen)],i->[ggen[i],fgen[i],1]));
    else
        todo := NewFIFO([[gone,fone,0]]);
    fi;
    if Length(ggen)=0 then # special case
        if compare(set,gone) then Add(result,fone); fi;
        return result;
    fi;
    seen := NewDictionary(ggen[1],false);
    forbidden := NewDictionary(Representative(f),false);
    for x in todo do
        if justone and KnowsDictionary(seen,x[1]) then
            AddDictionary(forbidden,x[2]);
            continue;
        fi;
        AddDictionary(seen,x[1]);
        if compare(set,x[1]) then Add(result,x[2]); if justone then break; fi; fi;
        if x[3]<n then
            for i in [1..Length(ggen)] do
                newfx := x[2]*fgen[i];
                newlen := x[3]+1;
                if Length(newfx)<newlen then continue; fi; # free reduction
                if ForAny([1..newlen],j->KnowsDictionary(forbidden,Subword(newfx,j,newlen))) then
                    continue;
                fi;
                Add(todo,[x[1]*ggen[i],newfx,newlen]);
            od;
            continue;
        fi;
    od;
    return result;
end);

InstallMethod(ShortGroupWordInSet, "(FR) for a group, an object and an int",
        [IsGroup,IsObject,IsObject],
        function(g,set,n)
    local f, s, sinv;
    s := GeneratorsOfGroup(g);
    f := FreeGroup(Length(s));
    sinv := Concatenation(List(s,x->[x^-1,x]));
    return SHORTWORDINSET@(f,GeneratorsOfMonoid(f),One(f),sinv,One(g),
            set,[GroupHomomorphismByImagesNC(f,g,GeneratorsOfGroup(f),s)],n);
end);

InstallMethod(ShortMonoidWordInSet, "(FR) for a monoid, an object and an int",
        [IsMonoid,IsObject,IsObject],
        function(g,set,n)
    local f, s;
    s := GeneratorsOfMonoid(g);
    f := FreeMonoid(Length(s));
    return SHORTWORDINSET@(f,GeneratorsOfMonoid(f),One(f),s,One(g),
            set,[FreeMonoidNatHomByGeneratorsNC(f,g)],n);
end);

InstallMethod(ShortSemigroupWordInSet, "(FR) for a semigroup, an object and an int",
        [IsSemigroup,IsObject,IsObject],
        function(g,set,n)
    local f, s;
    s := GeneratorsOfSemigroup(g);
    f := FreeSemigroup(Length(s));
    #!!! bug: GAP refuses free semigroups on 0 generators
    return SHORTWORDINSET@(f,GeneratorsOfSemigroup(f),fail,s,fail,
            set,[FreeSemigroupNatHomByGeneratorsNC(f,g)],n);
end);
#############################################################################

#############################################################################
##
#H Braid groups
##
InstallGlobalFunction(SurfaceBraidFpGroup, function(n,g,p)
  local G, R, a, b, z, s;
  a := List([1..g],i->CONCAT@("a",i));
  b := List([1..g],i->CONCAT@("b",i));
  s := List([1..n-1],i->CONCAT@("s",i));
  if p>0 then
    z := List([1..p-1],i->CONCAT@("z",i));
  else
    z := [];
  fi;
  G := FreeGroup(Concatenation(a,b,s,z));
  a := GeneratorsOfGroup(G){[1..g]};
  b := GeneratorsOfGroup(G){[g+1..2*g]};
  s := GeneratorsOfGroup(G){[2*g+1..2*g+n-1]};
  if p>0 then
    z := GeneratorsOfGroup(G){[2*g+n..2*g+n+p-2]};
  fi;
  R := [];
  Append(R,List(Filtered(Combinations([1..n-1],2),p->AbsInt(p[1]-p[2])>=2),p->Comm(s[p[1]],s[p[2]])));
  Append(R,List([1..n-2],i->s[i]*s[i+1]*s[i]/s[i+1]/s[i]/s[i+1]));
  Append(R,List(Cartesian([2..n-1],Concatenation(a,b)),p->Comm(s[p[1]],p[2])));
  if n>1 then
    Append(R,List(Concatenation(a,b),c->c*s[1]*c*s[1]/c/s[1]/c/s[1]));
    Append(R,List([1..g],i->a[i]*s[1]*b[i]/s[1]/a[i]/s[1]/b[i]/s[1]));
    Append(R,List(Cartesian([a,b],[a,b],Combinations([1..g],2)),p->Comm(p[1][p[3][2]],p[2][p[3][1]]^s[1])));
  fi;
  if p=0 then
    Add(R,Product([1..g],i->Comm(a[i],b[i]^-1),One(G))/Product([1..n-1],i->s[i],One(G))/Product([1..n-1],i->s[n-i],One(G)));
  else
    Append(R,List(Cartesian([2..n-1],z),p->Comm(s[p[1]],p[2])));
    Append(R,List(z,z->Comm(z,s[1]*z*s[1])));
    Append(R,List(Combinations(z,2),p->Comm(p[1]^s[1],p[2])));
    Append(R,List(Cartesian(z,Concatenation(a,b)),p->Comm(p[1]^s[1],p[2])));
  fi;
  return G/R;
end);

InstallGlobalFunction(PureSurfaceBraidFpGroup, function(n,g,p)
    local B, Bg, Sg, i;
    B := SurfaceBraidFpGroup(n,g,p);
    Bg := GeneratorsOfGroup(B);
    Sg := List([1..Length(Bg)],i->());
    for i in [1..n-1] do Sg[2*g+i] := (i,i+1); od;
    return Kernel(GroupHomomorphismByImages(B,SymmetricGroup(n),Bg,Sg));
end);

InstallGlobalFunction(CharneyBraidFpGroup, function(n)
    local B, Bg, Sg, i, f;
    B := SurfaceBraidFpGroup(n,0,1);
    Bg := GeneratorsOfGroup(B);
    Sg := List([1..n-1],i->(i,i+1));
    f := GroupHomomorphismByImages(B,SymmetricGroup(n),Bg,Sg);
    Sg := [];
    for i in SymmetricGroup(n) do if i<>() then
        Add(Sg,PreImagesRepresentativeNC(f,i));
    fi; od;
    return FpGroupPresentation(PresentationSubgroupMtc(B,Subgroup(B,Sg)));
end);

InstallGlobalFunction(ArtinRepresentation, function(n)
    local B, F, S;
    B := SurfaceBraidFpGroup(n,0,1);
    F := FreeGroup(n);
    S := GeneratorsOfGroup(F);
    return GroupHomomorphismByImages(B,AutomorphismGroup(F),
                   GeneratorsOfGroup(B),
                   List([1..n-1],i->MagmaEndomorphismByImagesNC(F,
                           Concatenation(S{[1..i-1]},
                                   [S[i+1],S[i]^S[i+1]],
                                   S{[i+2..n]}))));
end);

BindGlobal("ENDOISONE@", function(x)
    return ForAll(GeneratorsOfGroup(Source(x)),s->s=s^x);
end);

BindGlobal("ENDONORM@", function(x)
    if ENDOISONE@(x) then
        return 0;
    else
        return LogInt(Maximum(List(GeneratorsOfGroup(Source(x)),s->Length(s^x)))^4,2);
    fi;
end);
#############################################################################

#############################################################################
##
#M LowerCentralSeries etc. for algebras
##
InstallMethod(ProductIdeal, "for left ideal and right ideal",
        [IsAlgebra,IsAlgebra],
        function(i,j)
    local l, r, g;
    if HasLeftActingRingOfIdeal(i) then
        l := i;
    elif HasRightActingRingOfIdeal(i) then
        l := AsLeftIdeal(RightActingRingOfIdeal(i),i);
    else
        l := AsLeftIdeal(i,i);
    fi;
    if HasRightActingRingOfIdeal(j) then
        r := j;
    elif HasLeftActingRingOfIdeal(i) then
        r := AsRightIdeal(LeftActingRingOfIdeal(j),j);
    else
        r := AsRightIdeal(j,j);
    fi;
    g := [];
    for i in GeneratorsOfLeftIdeal(l) do
        for j in GeneratorsOfRightIdeal(r) do
            Add(g,i*j);
        od;
    od;
    return Ideal(LeftActingRingOfIdeal(l),g);
end);

InstallMethod(ProductBOIIdeal, "for two ideals over ring with invertible basis",
        [IsAlgebra,IsAlgebra],
        function(a,b)
    local g, r, i, j, k, c;
    if HasLeftActingRingOfIdeal(a) then
        r := LeftActingRingOfIdeal(a);
    elif HasLeftActingRingOfIdeal(b) then
        r := LeftActingRingOfIdeal(b);
    else
        r := a;
    fi;
    if HasGeneratorsOfIdeal(a) then
        a := GeneratorsOfIdeal(a);
    else
        a := GeneratorsOfAlgebra(a);
    fi;
    if HasGeneratorsOfIdeal(b) then
        b := GeneratorsOfIdeal(b);
    else
        b := GeneratorsOfAlgebra(b);
    fi;
    g := [];
    c := TwoSidedIdeal(r,g);
    for i in a do for j in b do
        k := i*j;
        if not k in c then
            Add(g,k);
            c := TwoSidedIdeal(r,g);
        fi;
    od; od;
    return c;
end);

BindGlobal("DIMENSIONSERIES@", function(A,n)
    local L, k, i;
    L := [AsTwoSidedIdeal(A,A)];
    i := AugmentationIdeal(A);
    k := i;
    while k<>L[Length(L)] and Length(L)<n do
        SetParent(k,A);
        Add(L,k);
        k := ProductBOIIdeal(k,i);
    od;
    return L;
end);

InstallMethod(DimensionSeries, "for an algebra with one",
        [IsAlgebra and HasAugmentationIdeal],
        A->DIMENSIONSERIES@(A,infinity));

InstallMethod(DimensionSeries, "for an algebra with one and a limit",
        [IsAlgebra and HasAugmentationIdeal,IsInt],
        DIMENSIONSERIES@);
#############################################################################

#############################################################################
##
#H Find incompressible elements
##
if false then
##!! this has nothing to do here
G := BinaryKneadingGroup(1/6);
S := [G.1,G.2,G.3];
pi := DecompositionOfFRElement(G);

EasyReduce := function(x)
  local e, i, verygeod, geod;
  e := ShallowCopy(ExtRepOfObj(x));
  verygeod := true;
  geod := true;
  for i in [2,4..Length(e)] do
    if e[i]=-1 then e[i] := 1; fi;
    if e[i] >= 2 or e[i] <= -2 then
      e[i] := RemInt(RemInt(e[i],2)+2,2);
      verygeod := false;
      if e[i-1] >= 2 then geod := false; fi;
    fi;
  od;
  return [ObjByExtRep(FamilyObj(x),e),geod,verygeod];
end;

MakeIncompressible := function(n)
  local inc, ginc, i;

  inc := [[One(G)],[G.1,G.2,G.3],Difference(Ball(G,2),Ball(G,1))];
  ginc := Ball(G,2);

  for i in [3..n] do
    inc[i+1] := Filtered(List(Cartesian(inc[i],[G.1,G.2,G.3]),p->p[1]*p[2]),function(g)
      local x;
      x := pi(g);
      return EasyReduce(g)[3] and x[1][1] in ginc and x[1][2] in ginc and EasyReduce(x[1][1])[2] and EasyReduce(x[1][2])[2];
    end);
    Append(ginc,inc[i+1]);
  od;
  return ginc;
end;
fi;
#############################################################################

#############################################################################
##
#F WreathProductPc
##
InstallOtherMethod( WreathProduct,
[IsPcGroup, IsPcGroup],
        function( G, H )
    return WreathProduct( G, H, RegularActionHomomorphism(H), Size(H));
end);

InstallOtherMethod( WreathProduct, true,
[IsPcGroup, IsPcGroup, IsMapping], 0,
        function( G, H, act )
    return WreathProduct( G, H, act,
                   Maximum( 1, LargestMovedPoint( Image( act ))));
end);

InstallOtherMethod( WreathProduct, true,
        [IsPcGroup, IsPcGroup, IsMapping, IsPosInt], 0,
function( G, H, act, l )
    local pcgsG, pcgsH, isoG, isoH, gensG, gensH, gensFG, gensFH, rels, W,
          i, j, k, m, n, a, b, F;

    pcgsG := Pcgs(G);
    pcgsH := Pcgs(H);
    n := Length( pcgsG );
    m := Length( pcgsH );

    F := FreeGroup(IsSyllableWordsFamily, m + n*l);
    rels := [];
    isoG := IsomorphismFpGroupByPcgs(pcgsG, "G");
    isoH := IsomorphismFpGroupByPcgs(pcgsH, "H");
    gensG := GeneratorsOfGroup(FreeGroupOfFpGroup(Range(isoG)));
    gensH := GeneratorsOfGroup(FreeGroupOfFpGroup(Range(isoH)));
    gensFG := List([1..l],i->GeneratorsOfGroup(F){m+(i-1)*n+[1..n]});
    gensFH := GeneratorsOfGroup(F){[1..m]};
    rels := List(RelatorsOfFpGroup(Range(isoH)),
                 w->MappedWord(w,gensH,gensFH));
    for i in [1..l] do
        Append(rels,List(RelatorsOfFpGroup(Range(isoG)),
                w->MappedWord(w,gensG,gensFG[i])));
    od;
    for i in [1..l] do
        for j in [i+1..l] do
            for a in gensFG[i] do
                for b in gensFG[j] do
                    Add(rels,Comm(a,b));
                od;
            od;
        od;
    od;
    for i in [1..l] do
        for a in [1..Length(pcgsH)] do
            b := pcgsH[a]^act;
            for j in [1..Length(gensFG[i])] do
                Add(rels,gensFG[i][j]^gensFH[a]/gensFG[i^b][j]);
            od;
        od;
    od;
    W := PcGroupFpGroup(F/rels);

    SetWreathProductInfo( W, rec(l := l, m := m, n := n,
        G := G, pcgsG := pcgsG, genG := GeneratorsOfGroup(G),
        H := H, pcgsH := pcgsH, genH := GeneratorsOfGroup(H),
        pcgsW := Pcgs(W),
        embeddings := []) );
    return W;
end );

InstallMethod(Embedding,"pc wreath product",
        [IsPcGroup and HasWreathProductInfo, IsPosInt],
        function(W,i)
    local info, FilledIn;

    FilledIn := function( exp, shift, len )
        local s;
        s := List([1..len], i->0);
        s{shift+[1..Length(exp)]} := exp;
        return s;
    end;

    info := WreathProductInfo(W);
    if not IsBound(info.embeddings[i]) then
        if i<=info.l then
            info.embeddings[i] := GroupHomomorphismByImagesNC(info.G,W,
                info.genG, List(info.genG, x->PcElementByExponents(info.pcgsW,
                    FilledIn(ExponentsOfPcElement(info.pcgsG,x),info.m+(i-1)*info.n,info.m+info.l*info.n))));
        elif i=info.l+1 then
            info.embeddings[i] := GroupHomomorphismByImagesNC(info.H,W,
                info.genH, List(info.genH, x->PcElementByExponents(info.pcgsW,
                    FilledIn(ExponentsOfPcElement(info.pcgsH,x),0,info.m+info.l*info.n))));
        else
            return fail;
        fi;
        SetIsInjective(info.embeddings[i],true);
    fi;
    return info.embeddings[i];
end);
#############################################################################

#############################################################################
# Dirichlet series
#############################################################################
BindGlobal("NEWDIRICHLETSERIES@", function(arg)
    return Objectify(NewType(DS_FAMILY,IsDirichletSeries),arg);
end);

InstallMethod(DirichletSeries, [], function()
    return NEWDIRICHLETSERIES@([],[],infinity);
end);

InstallMethod(DirichletSeries, [IsInt], function(n)
    return NEWDIRICHLETSERIES@([],[],n);
end);

InstallMethod(DirichletSeries, [IsList,IsList], function(ind,coeff)
    return NEWDIRICHLETSERIES@(ind,coeff,Maximum(ind));
end);

InstallMethod(DirichletSeries, [IsList,IsList,IsInt], function(ind,coeff,n)
    return NEWDIRICHLETSERIES@(ind,coeff,n);
end);
      
InstallMethod(DirichletSeries, [IsDirichletSeries,IsInt],
        function(s,n)
    local p;
    p := First([1..Length(s![1])+1],i->not IsBound(s![1][i]) or s![1][i] > n);
    return NEWDIRICHLETSERIES@(s![1]{[1..p-1]},s![2]{[1..p-1]},n);
end);

InstallMethod(ShrunkDirichletSeries, [IsDirichletSeries],
        function(s)
    local p, n;
    n := DegreeOfDirichletSeries(s);
    p := First([1..Length(s![1])+1],i->not IsBound(s![1][i]) or s![1][i] > n);
    return NEWDIRICHLETSERIES@(s![1]{[1..p-1]},s![2]{[1..p-1]},n);
end);

InstallMethod(String,[IsDirichletSeries],
        function(s)
    local n, first, v;
    v := "";
    for n in [1..Length(s![1])] do
        if s![2][n]<>0 then
            if s![2][n]>0 and v<>"" then
                Append(v,"+");
            fi;
            Append(v,String(s![2][n]));
            Append(v,"*"); Append(v,String(s![1][n])); Append(v,"^-s");
        fi;
    od;
    if v<>"" then Append(v,"+"); fi;
    Append(v,"o("); Append(v,String(s![3])); Append(v,"^-s)");
    return v;
end);

InstallMethod(ViewString, [IsDirichletSeries], String);

InstallMethod(DegreeOfDirichletSeries, [IsDirichletSeries],
        function(s)
    local i;
    i := Length(s![1]);
    while i>0 do
        if s![2][i]=0 then
            i := i-1;
        else
            return s![1][i];
        fi;
    od;
    return 0;
end);

InstallMethod(\+, [IsDirichletSeries,IsDirichletSeries],
        function(s1,s2)
    local i, p, maxdeg, coeff, val;
    
    maxdeg := Maximum(s1![3],s2![3]);
    
    coeff := ShallowCopy(s1![1]);
    val := ShallowCopy(s1![2]);
    
    for i in [1..Length(s2![1])] do
        p := PositionSorted(coeff,s2![1][i]);
        if p>Length(coeff) or coeff[p]<>s2![1][i] then
            Add(coeff,s2![1][i],p);
            Add(val,0,p);
        fi;
        val[p] := val[p] + s2![2][i];
    od;
    return NEWDIRICHLETSERIES@(coeff,val,maxdeg);
end);

InstallMethod(AdditiveInverseSameMutability, [IsDirichletSeries],
        function(s)
    return NEWDIRICHLETSERIES@(s![1],-s![2],s![3]);
end);

InstallMethod(Zero, [IsDirichletSeries],
        function(s)
    return NEWDIRICHLETSERIES@([],[],s![3]);
end);

InstallMethod(OneMutable, [IsDirichletSeries],
        function(s)
    return NEWDIRICHLETSERIES@([1],[1],s![3]);
end);

InstallMethod(\*, [IsDirichletSeries,IsScalar],
        function(s,x)
    return NEWDIRICHLETSERIES@(s![1],s![2]*x,s![3]);
end);

InstallMethod(\*, [IsScalar,IsDirichletSeries],
        function(x,s)
    return NEWDIRICHLETSERIES@(s![1],x*s![2],s![3]);
end);

InstallMethod(\*, [IsDirichletSeries,IsDirichletSeries],
        function(s1,s2)
    local coeff, val, i, j, degree, p, maxdeg;
    
    maxdeg := Maximum(s1![3],s2![3]);

    coeff := [];
    val := [];
 
    for i in [1..Length(s1![1])] do
        for j in [1..Length(s2![1])] do
            degree := s1![1][i]*s2![1][j];
            if degree > maxdeg then break; fi;
            p := PositionSorted(coeff,degree);
            if p>Length(coeff) or coeff[p]<>degree then
                Add(coeff,degree,p);
                Add(val,0,p);
            fi;
            val[p] := val[p] + s1![2][i]*s2![2][j];
        od;
    od;
    return NEWDIRICHLETSERIES@(coeff,val,maxdeg);
end);

InstallMethod(\=, [IsDirichletSeries,IsDirichletSeries],
        function(s1,s2)
    local i1, i2;
    if s1![3]<>s2![3] then return false; fi;
    i1 := 1;
    i2 := 1;
    while i1<=Length(s1![1]) or i2<=Length(s2![1]) do
        if i1<=Length(s1![1]) and (i2>Length(s2![1]) or s1![1][i1]<s2![1][i2]) then
            if s1![2][i1]<>0 then return false; fi;
            i1 := i1+1;
        elif i2<=Length(s2![1]) and (i1>Length(s1![1]) or s1![1][i1]>s2![1][i2]) then
            if s2![2][i2]<>0 then return false; fi;
            i2 := i2+1;
        else
            if s1![2][i1]<>s2![2][i2] then return false; fi;
            i1 := i1+1; i2 := i2+1;
        fi;
    od;
    return true;
end);

InstallMethod(\<, [IsDirichletSeries,IsDirichletSeries],
        function(s1,s2)
    local i1, i2;
    if s1![3]<s2![3] then return true; fi;
    if s1![3]>s2![3] then return false; fi;
    i1 := 1;
    i2 := 1;
    while i1<=Length(s1![1]) or i2<=Length(s2![1]) do
        if i1<=Length(s1![1]) and (i2>Length(s2![1]) or s1![1][i1]<s2![1][i2]) then
            if s1![2][i1]<0 then return true; fi;
            if s1![2][i1]>0 then return false; fi;
            i1 := i1+1;
        elif i2<=Length(s2![1]) and (i1>Length(s1![1]) or s1![1][i1]>s2![1][i2]) then
            if s2![2][i2]>0 then return true; fi;
            if s2![2][i2]<0 then return false; fi;
            i2 := i2+1;
        else
            if s1![2][i1]<s2![2][i2] then return true; fi;
            if s1![2][i1]>s2![2][i2] then return false; fi;
            i1 := i1+1; i2 := i2+1;
        fi;
    od;
    return false;
end);

InstallMethod(SpreadDirichletSeries, [IsDirichletSeries, IsInt],
        function(s,n)
    local p;
    p := First([1..Length(s![1])+1],i->not IsBound(s![1][i]) or s![1][i]^n > s![3]);
    return NEWDIRICHLETSERIES@(List(s![1]{[1..p-1]},i->i^n),s![2]{[1..p-1]},s![3]);
end);

InstallMethod(ShiftDirichletSeries, [IsDirichletSeries,IsInt],
        function(s,n)
    local p;
    p := First([1..Length(s![1])+1],i->not IsBound(s![1][i]) or s![1][i]*n > s![3]);
    return NEWDIRICHLETSERIES@(n*s![1]{[1..p-1]},s![2]{[1..p-1]},s![3]);
end);

InstallMethod(ZetaSeriesOfGroup, [IsGroup],
        function(G)
    local m;
    m := TransposedMat(CharacterDegrees(G));
    return DirichletSeries(m[1],m[2]);
end);

InstallMethod(ValueDirichletSeries, [IsDirichletSeries,IsRingElement],
        function(s,z)
    local v, i;
    v := Zero(z);
    for i in [1..Length(s![1])] do
        v := v+s![2][i]*s![1][i]^(-z);
    od;
    return v;
end);

InstallOtherMethod(Value, [IsDirichletSeries,IsRingElement],
        ValueDirichletSeries);
#############################################################################

#############################################################################
##
#F JenningsLieAlgebra
##
BindGlobal("LIEELEMENT@", function(A,l)
    return Objectify(A!.type,[A,l]);
end);

LIEEXTENDLCS@ := fail; # shut up warning

BindGlobal("LIECOMPUTEBASIS@", function(A,d)
    local i, l;

    if IsBound(A!.basis[d]) then return; fi;

    LIEEXTENDLCS@(A,d);
    A!.hom[d] := NaturalHomomorphismByNormalSubgroup(A!.lcs[d],A!.lcs[d+1]);
    A!.basis[d] := [];
    for i in IdentityMat(Length(AbelianInvariants(Range(A!.hom[d]))),A!.ring) do
        l := []; l[d] := i;
        Add(A!.basis[d],LIEELEMENT@(A,l));
    od;
    A!.transversal[d] := List(A!.pcp(Range(A!.hom[d])),
                              x->PreImagesRepresentativeNC(A!.hom[d],x));
end);

BindGlobal("JENNINGSSERIES@", function(G,p,d)
    local n, L, C, T;

    L := [G];
    for n in [2..d] do
        C := CommutatorSubgroup(G,L[n-1]);
        if p=0 then
            T := NaturalHomomorphismByNormalSubgroup(L[n-1],C);
            L[n] := PreImage(T,TorsionSubgroup(Range(T)));
        else
            L[n] := ClosureGroup(C,List(GeneratorsOfGroup(L[QuoInt(n+p-1,p)]), x->x^p));
        fi;
    od;
    return L;
end);

# a hack to shut up the GAP compiler warning
if not IsBound(PqEpimorphism) then PqEpimorphism := ReturnFail; fi;
if not IsBound(NqEpimorphismNilpotentQuotient) then NqEpimorphismNilpotentQuotient := ReturnFail; fi;

UnbindGlobal("LIEEXTENDLCS@");
BindGlobal("LIEEXTENDLCS@", function(A,d)
    local i;

    if d<=A!.degree then return; fi;

    if (IsBound(IsLpGroup) and IsLpGroup(A!.group)) or (IsFpGroup(A!.group) and Characteristic(A!.ring)=0) then
        A!.quo := NqEpimorphismNilpotentQuotient(A!.group,d);
        A!.pcp := Pcp;
        A!.exp := ExponentsByPcp;
    else
        A!.quo := PqEpimorphism(A!.group : ClassBound := d, Prime := Characteristic(A!.ring));
        A!.pcp := Pcgs;
        A!.exp := ExponentsOfPcElement;
    fi;
    A!.lcs := JENNINGSSERIES@(Range(A!.quo),Characteristic(A!.ring),d+1);
    A!.degree := d;
    for i in BoundPositions(A!.basis) do
        Unbind(A!.basis[i]);
        LIECOMPUTEBASIS@(A,i);
    od;
end);

if PqEpimorphism=ReturnFail then Unbind(PqEpimorphism); fi;
if NqEpimorphismNilpotentQuotient=ReturnFail then Unbind(NqEpimorphismNilpotentQuotient); fi;

InstallOtherMethod(JenningsLieAlgebra, "(FR) for a ring and an FP group",
        [IsRing,IsGroup],
        function(R,G)
    local C, A;

    C := NewFamily("LieFpElementsFamily",IsLieFpElementRep);
    A := Objectify(NewType(CollectionsFamily(C),
                 IsFpLieAlgebra and IsAttributeStoringRep),
                 rec(group := G,
                     family := C,
                     ring := R,
                     degree := 0,
                     lcs := [],
                     bracket := [],
                     pmap := [],
                     hom := [],
                     transversal := [],
                     basis := []));
    A!.type := NewType(A!.family,IsLieObject and IsLieFpElementRep);
    Grading(A);
    SetRepresentative(A,LIEELEMENT@(A,[]));
    SetLeftActingDomain(A,R);
    SetZero(C,Representative(A));
    return A;
end);

InstallMethod(GeneratorsOfAlgebra, "(FR) for an FP Lie algebra",
        [IsFpLieAlgebra],
        function(A)
    LIECOMPUTEBASIS@(A,1);
    return List(GeneratorsOfGroup(A!.group),x->LIEELEMENT@(A,[One(A!.ring)*A!.exp(A!.pcp(Range(A!.hom[1])),(x^A!.quo)^A!.hom[1])]));
end);

InstallMethod(Grading, "(FR) for a FP Lie algebra",
        [IsFpLieAlgebra],
        function(A)
    return rec(min_degree := 1, source := Integers, hom_components := function(d)
        LIECOMPUTEBASIS@(A,d);
        return VectorSpace(A!.ring,A!.basis[d],Representative(A));
    end);
end);

InstallMethod(ViewString, "(FR) for a FP Lie algebra",
        [IsFpLieAlgebra],
        function(A)
    local s;
    s := Concatenation("<FP Lie algebra over ",String(A!.ring));
    if A!.degree>0 then
        if IsTrivial(A!.lcs[A!.degree]) then
            Append(s,", of");
        else
            Append(s,", computed up to");
        fi;
        APPEND@(s," degree ",A!.degree);
    fi;
    Append(s,">");
    return s;
end);
INSTALLPRINTERS@(IsFpLieAlgebra);

InstallOtherMethod(ZeroOp, "(FR) for a FP Lie algebra",
        [IsFpLieAlgebra],
        function(A)
    return Representative(A);
end);

InstallMethod(ZeroOp, "(FR) for a FP Lie algebra element",
        [IsLieObject and IsLieFpElementRep],
        function(X)
    return LIEELEMENT@(X![1],[]);
end);

InstallMethod(EQ, "(FR) for FP Lie algebra elements",
        IsIdenticalObj,
        [IsLieObject and IsLieFpElementRep,IsLieObject and IsLieFpElementRep],
        function(X,Y)
    local n;
    for n in Union(BoundPositions(X![2]),BoundPositions(Y![2])) do
        if not IsBound(X![2][n]) and not ForAll(Y![2][n],IsZero) then
            return false;
        fi;
        if not IsBound(Y![2][n]) and not ForAll(X![2][n],IsZero) then
            return false;
        fi;
        if X![2][n]<>Y![2][n] then
            return false;
        fi;
    od;
    return true;
end);

InstallMethod(LT, "(FR) for FP Lie algebra elements",
        IsIdenticalObj,
        [IsLieObject and IsLieFpElementRep,IsLieObject and IsLieFpElementRep],
        function(X,Y)
    local n;
    for n in Union(BoundPositions(X![2]),BoundPositions(Y![2])) do
        if not IsBound(X![2][n]) and not ForAll(Y![2][n],IsZero) then
            return Y![2][n]>0;
        fi;
        if not IsBound(Y![2][n]) and not ForAll(X![2][n],IsZero) then
            return X![2][n]<0;
        fi;
        if X![2][n]<>Y![2][n] then
            return X![2][n]<Y![2][n];
        fi;
    od;
    return false;
end);

InstallMethod(\+, "(FR) for FP Lie algebra elements",
        IsIdenticalObj,
        [IsLieObject and IsLieFpElementRep,IsLieObject and IsLieFpElementRep],
        function(X,Y)
    local m, n;
    m := [];
    for n in Union(BoundPositions(X![2]),BoundPositions(Y![2])) do
        if not IsBound(X![2][n]) then
            m[n] := Y![2][n];
        elif not IsBound(Y![2][n]) then
            m[n] := X![2][n];
        else
            m[n] := X![2][n]+Y![2][n];
        fi;
    od;
    return LIEELEMENT@(X![1],m);
end);

InstallMethod(\-, "(FR) for FP Lie algebra elements",
        IsIdenticalObj,
        [IsLieObject and IsLieFpElementRep,IsLieObject and IsLieFpElementRep],
        function(X,Y)
    local m, n;
    m := [];
    for n in Union(BoundPositions(X![2]),BoundPositions(Y![2])) do
        if not IsBound(X![2][n]) then
            m[n] := -Y![2][n];
        elif not IsBound(Y![2][n]) then
            m[n] := X![2][n];
        else
            m[n] := X![2][n]-Y![2][n];
        fi;
    od;
    return LIEELEMENT@(X![1],m);
end);

InstallMethod(AdditiveInverseSameMutability, "(FR) for an FP Lie algebra element",
        [IsLieObject and IsLieFpElementRep],
        function(X)
    local m, n;
    m := [];
    for n in BoundPositions(X![2]) do
        m[n] := -X![2][n];
    od;
    return LIEELEMENT@(X![1],m);
end);

InstallMethod(\*, "(FR) for FP Lie algebra elements",
        IsIdenticalObj,
        [IsLieObject and IsLieFpElementRep,IsLieObject and IsLieFpElementRep],
        function(X,Y)
    local i, j, ii, jj, m, A;
    m := [];
    A := X![1];
    for i in BoundPositions(X![2]) do
        if not IsBound(A!.bracket[i]) then A!.bracket[i] := []; fi;
        for j in BoundPositions(Y![2]) do
            if not IsBound(A!.bracket[i][j]) then A!.bracket[i][j] := []; fi;
            LIECOMPUTEBASIS@(A,i+j);
            if not IsBound(m[i+j]) then
                m[i+j] := List(A!.basis[i+j],i->Zero(A!.ring));
            fi;
            for ii in [1..Length(A!.basis[i])] do
                if not IsBound(A!.bracket[i][j][ii]) then
                    A!.bracket[i][j][ii] := [];
                fi;
                for jj in [1..Length(A!.basis[j])] do
                    if not IsBound(A!.bracket[i][j][ii][jj]) then
                        A!.bracket[i][j][ii][jj] := A!.exp(A!.pcp(Range(A!.hom[i+j])),Comm(A!.transversal[i][ii],A!.transversal[j][jj])^A!.hom[i+j]);
                    fi;
                    m[i+j] := m[i+j] + X![2][i][ii]*Y![2][j][jj]*A!.bracket[i][j][ii][jj];
                od;
            od;
        od;
    od;
    for i in BoundPositions(m) do
        if ForAll(m[i],IsZero) then Unbind(m[i]); fi;
    od;
    return LIEELEMENT@(X![1],m);
end);

InstallMethod(\*, "(FR) for a scalar and an FP Lie algebra element",
        [IsScalar,IsLieObject and IsLieFpElementRep],
        function(X,Y)
    return LIEELEMENT@(Y![1],X*Y![2]);
end);

InstallMethod(\*, "(FR) for an FP Lie algebra element and a scalar",
        [IsLieObject and IsLieFpElementRep,IsScalar],
        function(X,Y)
    return LIEELEMENT@(X![1],X![2]*Y);
end);

BindGlobal("PTHPOWER@", function(X,A,p,s)
    local m, i, ii, j, t;
    m := [];
    for i in BoundPositions(X![2]) do
        LIECOMPUTEBASIS@(A,p*i);
        if not IsBound(m[p*i]) then
            m[p*i] := List(A!.basis[p*i],i->Zero(A!.ring));
        fi;
        if not IsBound(A!.pmap[i]) then A!.pmap[i] := []; fi;
        for ii in [1..Length(A!.basis[i])] do
            if not IsBound(A!.pmap[i][ii]) then
                A!.pmap[i][ii] := A!.exp(A!.pcp(Range(A!.hom[p*i])),(A!.transversal[i][ii]^p)^A!.hom[p*i]);
            fi;
            m[p*i] := m[p*i] + X![2][i][ii]^p*A!.pmap[i][ii];
        od;
    od;
    m := LIEELEMENT@(X![1],m);
    for i in BoundPositions(X![2]) do
        for ii in [1..Length(A!.basis[i])] do
            if IsBound(X![2][i]) and not IsZero(X![2][i][ii]) then
                t := X![2][i][ii]*A!.basis[i][ii];
                X := X - t;
                for j in [1..p-1] do m := m + s[j]([t,X]); od;
            fi;
        od;
    od;
    return m;
end);

InstallMethod(\^, "(FR) for an FR Lie algebra element and a p-power",
        [IsLieObject and IsLieFpElementRep,IsPosInt],
        function(X,Y)
    local p, n;
    p := Characteristic(X![1]!.ring);
    if p=0 or p^LogInt(Y,p)<>p then
        Error(Y," must be a power of the characteristic");
    fi;
    for n in [1..LogInt(Y,p)] do
        X := PTHPOWER@(X,X![1],p,PowerS(X![1]));
    od;
    return X;
end);

InstallMethod(Degree, "(FR) for a FR Lie algebra element",
        [IsLieObject and IsLieFpElementRep],
         function(X)
    local n;
    for n in [1..Length(X![2])] do
        if IsBound(X![2][n]) and not ForAll(X![2][n],IsZero) then return n; fi;
    od;
    return infinity;
end);

InstallMethod(ViewString, "(FR) for a FP Lie algebra element",
        [IsLieObject and IsLieFpElementRep],
        function(X)
    local n, s;
    s := "<";
    n := Degree(X);
    if n=infinity then
        Append(s,"zero");
    else
        if ForAll([n+1..Length(X![2])],i->not IsBound(X![2][i]) or ForAll(X![2][i],IsZero)) then
            Append(s,"homogeneous ");
        fi;
        APPEND@(s,"degree-",n);
    fi;
    Append(s," Lie element>");
    return s;
end);
INSTALLPRINTERS@(IsLieObject and IsLieFpElementRep);

BindGlobal("LIE2VECTOR@", function(dims,dim,x)
    local i, v;
    v := List([1..dim],i->Zero(x![1]!.ring));
    for i in BoundPositions(x![2]) do
        if not IsBound(dims[i]) then
            if not ForAll(x![2][i],IsZero) then
                return fail;
            fi;
        else
            v{dims[i]} := x![2][i];
        fi;
    od;
    ConvertToVectorRep(v);
    MakeImmutable(v);
    return v;
end);

InstallHandlingByNiceBasis("IsLieFpElementSpace", rec(
        detect := function(R,l,V,z)
    if l=[] then
        return IsLieObject(z) and IsLieFpElementRep(z) and z![1]!.ring=R;
    else
        return ForAll(l,x->IsLieObject(x) and IsLieFpElementRep(x) and x![1]!.ring=R);
    fi;
end,
  NiceFreeLeftModuleInfo := function(V)
    local b, dim, dims, i, j, k, R;
    b := GeneratorsOfLeftModule(V);
    R := Zero(V)![1];
    dim := 0;
    dims := [];
    for i in b do
        for j in BoundPositions(i![2]) do
            if ForAll(i![2][j],IsZero) then continue; fi;
            dims[j] := true;
        od;
    od;
    for i in BoundPositions(dims) do
        dims[i] := [dim+1..dim+Length(R!.basis[i])];
        dim := dim+Length(R!.basis[i]);
    od;
    b := List(b,x->LIE2VECTOR@(dims,dim,x));
    return rec(dims := dims,
               dim := dim,
               ring := R,
               space := VectorSpace(LeftActingDomain(V),b,LIE2VECTOR@(dims,dim,Zero(V))));
end,
  NiceVector := function(V,v)
    local info, x;
    info := NiceFreeLeftModuleInfo(V);
    x := LIE2VECTOR@(info.dims,info.dim,v);
    if x=fail or not x in info.space then
        return fail;
    else
        return x;
    fi;
end,
  UglyVector := function(V,v)
    local info, i, l;
    info := NiceFreeLeftModuleInfo(V);
    if not v in info.space then return fail; fi;
    l := [];
    for i in BoundPositions(info.dims) do
        l[i] := v{info.dims[i]};
    od;
    return LIEELEMENT@(info.ring,l);
end));
#############################################################################

# solve linear system mod N.
# returns x such that x*mat = vec mod N, or "fail".
# note that x could be rational, non-integer.
InstallMethod(SolutionMatModN, [IsMatrix,IsVector,IsPosInt],
        function(mat,vec,N)
    local sol, i, p, s0, M, row;
    
    sol := List(mat,row->0);
    if N=1 then return sol; fi; # bug with FactorsInt containing 1
    
    mat := SmithNormalFormIntegerMatTransforms(mat);
    vec := vec*mat.coltrans;
    row := mat.rowtrans;
    mat := mat.normal;
    for i in [1..Minimum(Length(mat),Length(mat[1]))] do
        p := AbsoluteValue(mat[i][i]);
        if p>1 then
            row[i] := row[i] / p;
            mat[i][i] := mat[i][i] / p;
        fi;
    od;

    M := 1;
    for p in FactorsInt(N) do
        s0 := SolutionMat(mat*Z(p),vec*Z(p));
        if s0=fail then return fail; fi;
        s0 := List(s0,IntFFE);
        vec := (vec - s0*mat)/p;
        sol := sol + M*s0;
        M := M*p;
    od;
    
    return sol*row;
end);

BindGlobal("FRAC@", function(x)
    x := x-Int(x);
    if x>=0 then return x; else return x+1; fi;
end);

# solve rational linear equation in Q/Z.
# for now, only treat integer matrix.

# overlaps with SolveHomEquationsModZ in package Cryst.
InstallMethod(SolutionMatMod1, [IsMatrix, IsVector],
        function(mat,vec)
    local sol, N, d, i, snf, row;
    
    # non-optimal: should store the Smith normal form as an attribute
    snf := SmithNormalFormIntegerMatTransforms(mat);
    vec := vec*snf.coltrans;
    row := ShallowCopy(snf.rowtrans);
    sol := [];

    for i in [1..Length(vec)] do
        # diagonal term
        if i<=Length(snf.normal) then
            d := snf.normal[i][i];
        else
            d := 0;
        fi;
        
        if d=0 then
            if FRAC@(vec[i])<>0 then return fail; fi; # no solution
            
            if i<=Length(snf.normal) then
                Add(sol,[0]); # infinite set of solutions, this is one
            fi;
        else
            Add(sol,(vec[i]+[0..d-1])/d);
        fi;
    od;
    return Set(Cartesian(sol)*row,x->List(x,FRAC@));
end);

# argument of a cyclotomic number, assumed to be a multiple of a root of unity.
# returns a rational number in [0,1).
InstallMethod(CyclotomicByArgument, [IsRat], function(q)
    return E(DenominatorRat(q))^NumeratorRat(q);
end);

InstallMethod(ArgumentOfCyclotomic, [IsCyc], function(z)
    local q;
    q := DescriptionOfRootOfUnity(z);
    q := q[2]/q[1];
    q := q-Int(q);
    if z<>CyclotomicByArgument(q) then
        TryNextMethod(); # do floating-point approximations, probably
    fi;
    return q;
end);

################################################################
InstallMethod(IrreducibleRepresentations@, [IsGroup], function(G)
    local reps, r;
    reps := IrreducibleRepresentations(G);
    for r in reps do SetIsLinearRepresentation(r,true); od;
    return reps;
end);

InstallMethod(ProjectiveRepresentationByFunction, [IsGroup,IsGroup,IsFunction],
        function(g,h,f)
    local r;
    r := MappingByFunction(g,h,f);
    SetIsProjectiveRepresentation(r,true);
    return r;
end);

InstallMethod(LinearRepresentationByImages, [IsGroup,IsGroup,IsList,IsList],
        function(g,h,src,img)
    local r;
    r := GroupHomomorphismByImages(g,h,src,img);
    SetIsLinearRepresentation(r,true);
    return r;
end);

InstallMethod(DegreeOfProjectiveRepresentation, [IsProjectiveRepresentation],
        function(rep)
    return Length(Image(rep,One(Source(rep))));
end);
        
InstallMethod(Degree, [IsProjectiveRepresentation],
        DegreeOfProjectiveRepresentation);

# the projective representation extending the linear representation "rep".
# returns a group homomorphism (which is not really a homomorphism! only
# up to scalars) with source "group".
InstallMethod(ProjectiveExtension, [IsLinearRepresentation, IsGroup],
        function(rep, group)
    local n, rank, transversal, mat, res, t, g, P, m, a, b;
     
    n := Source(rep);
    
    if n=group then # special cases, which GAP doesn't handle well
        return rep;
    elif IsTrivial(n) then
        return LinearRepresentationByImages(group,Range(rep),GeneratorsOfGroup(group),List(GeneratorsOfGroup(group),x->One(Range(rep))));
    fi;
    
    rank := Length(Image(rep,One(n)));
 
    res := [];
    transversal := List(RightTransversal(group,n),x->CanonicalRightCosetElement(n,x));
    for t in transversal do
        P := [];
        for g in GeneratorsOfGroup(n) do
            a := Image(rep, g);
            b := Image(rep, g^t);
            for m in Basis(MatrixAlgebra(Rationals,rank)) do
                Add(P, Concatenation(TransposedMat(a)*m - m*TransposedMat(b)));
            od;
        od;
        mat := NullspaceMat(TransposedMat(P));
        if Length(mat)<>1 then
            Error("Solution is not 1-dimensional!");
        fi;
        mat := List([0..rank-1], i->mat[1]{rank*i+[1..rank]});
        if AbsoluteValue(Norm(DeterminantMat(mat)))<>1 then
            Print("# warning: got a matrix with determinant not of norm 1\n");
            mat := mat / DeterminantMat(mat)^(1/rank);
        fi;
        Add(res, mat);
    od;
    return ProjectiveRepresentationByFunction(group,Range(rep),function(x)
        local t;
        t := CanonicalRightCosetElement(n,x);
        return Image(rep,x/t)*res[Position(transversal,t)];
    end);
end);    

InstallMethod(ProjectiveExtension, [IsProjectiveRepresentation, IsGroup],
        function(rep, group)
    # one easy case
    if Source(rep)=group then
        return rep;
    fi;
    
    # don't bother
    TryNextMethod();
end);

InstallMethod(ProjectiveQuotient, [IsProjectiveRepresentation,IsGroupHomomorphism],
        function(rep, epi)
    return ProjectiveRepresentationByFunction(Image(epi),Range(rep),
                   x->CanonicalRightCosetElement(Kernel(epi),PreImagesRepresentativeNC(epi,x))^rep);
end);
        
InstallMethod(CoboundaryMatrix, [IsGroup], function(G)
    local mat, i, j, k, m, elements;
    
    elements := AsSortedList(G);
    mat := [];
    for i in [1..Length(elements)] do
        m := NullMat(Length(elements),Length(elements));
        for j in [1..Length(elements)] do
            m[i][j] := m[i][j] + 1;
            m[j][i] := m[j][i] + 1;
            k := Position(elements,elements[j]^-1*elements[i]);
            m[j][k] := m[j][k] - 1;
        od;
        Add(mat,Concatenation(m));
    od;
    return mat;
end);

InstallMethod(EpimorphismSchurCover@, [IsGroup], EpimorphismSchurCover);

InstallMethod(EpimorphismSchurCover@, [IsPcGroup],
        function(G)
    local c;
    c := EpimorphismSchurCover(G);
    c := InverseGeneralMapping(IsomorphismPcGroup(Source(c)))*c;
    IsGroupHomomorphism(c);
    return c;
end);

InstallMethod(EpimorphismSchurCover@, [IsPermGroup],
        function(G)
    local c;
    c := EpimorphismSchurCover(G);
    return InverseGeneralMapping(IsomorphismPermGroup(Source(c)))*c;
end);

InstallMethod(TensorProductOp, [IsLinearRepresentation,IsLinearRepresentation],
        function(r1,r2)
    local range, g, gens, img;
    
    g := Source(r1);
    gens := GeneratorsOfGroup(g);
    img := List(gens,x->KroneckerProduct(x^r1,x^r2));
    range := GroupByGenerators(img,KroneckerProduct(One(Range(r1)),One(Range(r2))));
    return LinearRepresentationByImages(g,range,gens,img);
end);

InstallMethod(TensorProductOp, [IsProjectiveRepresentation,IsProjectiveRepresentation],
        function(r1,r2)
    local range, g, img;
    
    g := Source(r1);
    img := List(AsSortedList(g),x->KroneckerProduct(x^r1,x^r2));
    
    return ProjectiveRepresentationByFunction(g,Group(img),x->KroneckerProduct(x^r1,x^r2));
end);

# compares two cocycles QxQ->Q/Z, and says if they differ by a coboundary.
InstallMethod(AreCohomologous, [IsList,IsList,IsGroup],
        function(c1,c2,Q)
    local elements, diff, i, j, k, m, denom;
    
    if c1=c2 then return true; fi; # speedup
    
    return SolutionMatMod1(CoboundaryMatrix(Q),c1-c2)<>fail;
end);

