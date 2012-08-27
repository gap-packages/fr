# todo:
# 1- refine meshes, esp. at pts of large degree
# 2- determine formula for edge shrinking wrt degrees
# 3- test

tri := TRIVIALSPIDER@FR([P1infinity,P1one,P1zero]);
g := FreeGroup(3);
IMGMARKING@FR(tri,g);
#permrep := GroupHomomorphismByImages(g,SymmetricGroup(5),GeneratorsOfGroup(g),[(1,2,3,4,5),(1,2),((1,2,3,4,5)*(1,2))^-1]);
permrep := GroupHomomorphismByImages(g,SymmetricGroup(3),GeneratorsOfGroup(g),[(1,2,3),(1,2),(2,3)]);


liftspider := function(spider,monodromy)
    # create a new triangulation lifted by the given monodromy rep'n.
    # the positions in the triangulation are the same as in the original
    # spider; i.e., geometrically one takes deg(monodromy) copies of the
    # original sphere, cuts it along the minimal spanning tree, and glues
    # the sheets to each other along the monodromy rep'n.
    # vertices of the result acquire an extra field, "degree", which is the
    # local degree of the map.

    local d, e, f, v, i, j, edgeperm, reverse,
          edges, faces, lift;
    
    Assert(0,Source(spider!.marking)=Source(monodromy));
    
    d := Maximum(LargestMovedPoint(Range(monodromy)),1);
    Assert(0,IsTransitive(Image(monodromy),[1..d]));
    Assert(0,SPIDERRELATOR@FR(spider)^monodromy=());

    # make d copies of the faces and edges, renumber their indices
    faces := [];
    edges := List([1..d],i->StructuralCopy(spider!.cut!.e));
    for i in [1..d] do
        faces[i] := [];
        for e in edges[i] do
            e.index := e.index + (i-1)*Length(spider!.cut!.e);
            j := e.left.index;
            if j<=Length(spider!.cut!.f) and not IsBound(faces[i][j]) then
                faces[i][j] := e.left;
                e.left.index := j + (i-1)*Length(spider!.cut!.f);
            fi;
        od;
    od;
    
    lift := rec(v := [],
                e := Concatenation(edges),
                f := Concatenation(faces));
    
    # reattach the edges according to the monodromy action
    edgeperm := List(edges[1],e->PreImagesRepresentative(spider!.marking,e.gpelement)^monodromy);
    reverse := List(edges[1],e->e.reverse.index);
    for i in [1..d] do
        for j in [1..Length(edges[i])] do
            edges[i][j].reverse := edges[i^edgeperm[j]][reverse[j]];
            edges[i][j].right := edges[i][j].reverse.left;
        od;
    od;
    
    # create new vertices, attach the edges to them
    for e in lift.e do
        if e.from.type='v' then # old vertex, replace it simply
            v := rec(type := 'w', # temporary other letter
                     n := [],
                     pos := e.from.pos, # keep old positions for a moment
                     degree := 1/Length(e.from.n),
                     index := Length(lift.v)+1,
                     operations := edges[1][1].operations);
            Add(lift.v,v);
            if IsBound(spider!.cut!.v[i].fake) then
                v.fake := true;
            fi;
            repeat
                e.from := v;
                Add(v.n,e);
                i := POSITIONID@FR(e.left.n,e)-1;
                if i=0 then i := Length(e.left.n); fi;
                e := e.left.n[i].reverse;
            until IsIdenticalObj(e,v.n[1]);
            v.degree := v.degree * Length(v.n);
        fi;
    od;
    
    for v in lift.v do
        e.type := 'v';
    od;
    
    # correct to pointers
    for e in lift.e do
        e.to := e.reverse.from;
    od;
    
    return Objectify(TYPE_TRIANGULATION,lift);
end;

refinetriangulation := function(triangulation,maxlen)
    # refines triangulation by adding circumcenters, until
    # the length of every edge (say connecting v to w) is at most
    # maxlen^Maximum(v.degree*w.degree)
    local idle, e, f, maxdegree, len, mult, p;
    
    maxdegree := Maximum(List(triangulation!.v,v->v.degree));
    mult := Int(Log(7./maxlen)/Log(1.5));
    repeat
        len := List([1..maxdegree],i->(maxlen*1.5^mult)^i);
        idle := true;
        for e in triangulation!.e do
            if P1Distance(e.from.pos,e.to.pos) > len[Maximum(e.from.degree,e.to.degree)] then
                idle := false;
                f := e.left;
                if not IsBound(f.radius) then
                    p := CallFuncList(P1Circumcentre,List(f.n,e->e.from.pos));
                    f.centre := p[1];
                    f.radius := p[2];
                fi;
                ADDTOTRIANGULATION@FR(triangulation,f,f.centre);
                triangulation!.v[Length(triangulation!.v)].degree := 1;
            fi;
        od;
        if idle then mult := mult-1; fi;
    until idle and mult<0;
    
    # that's a waste of time, we won't use these values
    for e in triangulation!.e do
        if not IsBound(e.pos) then
            e.pos := P1Barycentre(e.from.pos,e.to.pos);
            e.map := EDGEMAP@FR(e);
        fi;
    od;
    for f in triangulation!.f do
        if not IsBound(f.pos) then
            f.pos := P1Barycentre(List(f.n,e->e.to.pos));
        fi;
    od;
end;

layouttriangulation := function(triangulation)
    # run matlab code to optimize point placement
    local i, e, f, v, m, max, infty, tmpin, tmpout, dir, file;
    
    dir := DirectoryTemporary();
    tmpin := Filename(dir,"triangulation");
    tmpout := Filename(dir,"positions");
    file := OutputTextFile(tmpin,false);
    
    max := 0;
    for v in triangulation!.v do
        m := Length(v.n);
        if m>max then
            infty := v.index;
            max := m;
        fi;
    od;
    PrintTo(file,"VERTICES ",Length(triangulation!.v),"\n",infty,"\n");
    
    PrintTo(file,"EDGES ",Length(triangulation!.e),"\n");
    for e in triangulation!.e do
        PrintTo(file,e.from.index," ",e.to.index," ",P1Distance(e.from.pos,e.to.pos),"\n");
        #!!! distance should be something like (L/2)^(1/e.from.degree)+(L/2)^(1/e.to.degree)
    od;
    
    PrintTo(file,"FACES ",Length(triangulation!.f),"\n");
    for f in triangulation!.f do
        for m in [1..3] do
            PrintTo(file,f.n[m].from.index);
            if m=3 then PrintTo(file,"\n"); else PrintTo(file," "); fi;
        od;
    od;
    CloseStream(file);
    Process(DirectoryCurrent(),"/usr/local/bin/matlab",
            InputTextNone(),OutputTextUser(),
            ["-nosplash","-nojvm",
             "-r",Concatenation("layout('",tmpin,"','",tmpout,"');")]);
    file := InputTextFile(tmpout);
    m := ReadAll(file);
    CloseStream(file);
    m := EvalString(m);
    Remove(m);
    
    v := List(m,P1Sphere);
    m := NORMALIZINGMAP@FR(v,fail);
    v := List(v,v->P1Image(m,v));
    
    for i in [1..Length(v)] do
        triangulation!.v[i].pos := v[i];
    od;
    for e in triangulation!.e do
        e.pos := P1Midpoint(e.from.pos,e.to.pos);
    od;
    for f in triangulation!.f do
        f.pos := P1Barycentre(List(f.n,x->x.from.pos));
        i := CallFuncList(P1Circumcentre,List(f.n,e->e.from.pos));
        f.centre := i[1];
        f.radius := i[2];
    od;
end;

t := liftspider(tri,permrep);
refinetriangulation(t,0.3);
layouttriangulation(t);
#Draw(t);