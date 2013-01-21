EpimorphismGermGroup := function(G,n)
    local P, pcpP, pcpG, Q, l, i, data, act;
    P := PermGroup(G,n);
    i := IsomorphismPcpGroup(P);
    if i=fail then return fail; fi;
    data := GermData(G);
    l := Length(Alphabet(G))^n;
    Q := WreathProduct(data.group,Range(i),InverseGeneralMapping(i),l);
    pcpP := Pcp(Range(i));
    pcpG := Pcp(data.group);
    act := function(g)
        local e, d, x;
        d := Decomposition(g,n);
        e := ExponentsByPcp(pcpP,PermList(d[2])^i);
        for x in d[1] do
            Append(e,ExponentsByPcp(pcpG,GermValue(x,data)));
        od;
        return PcpElementByExponents(Collector(Q),e);
    end;
    return GroupHomomorphismByFunction(G,Q,act);
end;

# stolen from pkg/polycyclic/gap/pcpgrp/wreath.gi 1.1
InstallOtherMethod( WreathProduct, true, 
[IsPcpGroup, IsPcpGroup, IsMapping, IsInt], 0,
function( G, H, act, l )
    local pcpG, relG, pcpH, relH, n, m, coll, i, k, c, e, o, j, W, a, h,
          ShiftedObject;

    ShiftedObject := function( exp, c )
        local obj, i;
        obj := [];
        for i in [1..Length(exp)] do
            if exp[i] <> 0 then Append( obj, [c+i, exp[i]] ); fi;
        od;
        return obj;
    end;

    pcpG := Pcp(G);
    relG := RelativeOrdersOfPcp( pcpG );
    pcpH := Pcp(H);
    relH := RelativeOrdersOfPcp( pcpH );
    n := Length( pcpG );
    m := Length( pcpH );

    coll := FromTheLeftCollector( m + n*l );

    # relations of G
    for i in [1..n] do

        if relG[i] > 0 then
            e := ExponentsByPcp( pcpG, pcpG[i]^relG[i] );
            for k in [1..l] do
                c := m + (k-1)*n;
                o := ShiftedObject( e, c );
                SetRelativeOrder( coll, c+i, relG[i] );
                SetPower( coll, c+i, o );
            od;
        fi;

        for j in [1..i-1] do
            e := ExponentsByPcp( pcpG, pcpG[i]^pcpG[j] );
            for k in [1..l] do
                c := m + (k-1)*n;
                o := ShiftedObject( e, c );
                SetConjugate( coll, c+i, c+j, o );
            od;

            e := ExponentsByPcp( pcpG, pcpG[i]^(pcpG[j]^-1) );
            for k in [1..l] do
                c := m + (k-1)*n;
                o := ShiftedObject( e, c );
                SetConjugate( coll, c+i, -(c+j), o );
            od;
        od;
    od;

    # relations of H
    for i in [1..m] do

        if relH[i] > 0 then
            e := ExponentsByPcp( pcpH, pcpH[i]^relH[i] );
            o := ShiftedObject( e, 0 );
            SetRelativeOrder( coll, i, relH[i] );
            SetPower( coll, i, o );
        fi;

        for j in [1..i-1] do

            e := ExponentsByPcp( pcpH, pcpH[i]^pcpH[j] );
            o := ShiftedObject( e, 0 );
            SetConjugate( coll, i, j, o );

            e := ExponentsByPcp( pcpH, pcpH[i]^(pcpH[j]^-1) );
            o := ShiftedObject( e, 0 );
            SetConjugate( coll, i, -j, o );

        od;
    od;

    # action of H
    for j in [1..m] do
        a := Image( act, pcpH[i] );
        for k in [1..l] do
            h := k^a;
            for i in [1..n] do
                o := [m + (h-1)*n + i, 1];
                SetConjugate( coll, m + (k-1)*n + i, j, o );
            od;
            h := k^(a^-1);
            for i in [1..n] do
                o := [m + (h-1)*n + i, 1];
                SetConjugate( coll, m + (k-1)*n + i, -j, o );
            od;
        od;
    od;

    UpdatePolycyclicCollector( coll );
    W := PcpGroupByCollectorNC( coll );
    return W;
end );

