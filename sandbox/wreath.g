DeclareOperation("WreathProductOp", [IsGroup, IsGroup, IsGroupHomomorphism, IsListOrCollection]);
MakeReadWriteGlobal("WreathProduct");

WreathProduct := function(arg)
    if Length(arg)=2 and IsPermGroup(arg[2]) then
        return WreathProductImprimitiveAction(arg[1],arg[2]);
    elif Length(arg)=2 and IsExternalSet(arg[2]) then
        return WreathProductOp(arg[1],ActingDomain(arg[2]),ActionHomomorphism(arg[2]),arg[2]);
    elif Length(arg)=2 and IsGroup(arg[2]) then
        return WreathProductOp(arg[1],arg[2],RegularActionHomomorphism(arg[2]),[1..Size(arg[2])]);
    elif Length(arg)=3 then
        return WreathProductOp(arg[1],arg[2],arg[3],[1..LargestMovedPoint(Image(arg[3]))]);
    elif Length(arg)=4 and IsInt(arg[4]) then
        return WreathProductOp(arg[1],arg[2],arg[3],[1..arg[4]]);
    elif Length(arg)=4 and IsList(arg[4]) then
        return WreathProductOp(arg[1],arg[2],arg[3],arg[4]);
    else
        Error("Invalid arguments to WreathProduct");
    fi;
end;

InstallMethod(WreathProductOp, "for pc groups",
#XXX: Write the more efficient p-group version (this is just the
#     soluble version, but is already ridiculously faster than the
#     perm version for Sylow subgroups of Sym(p^n).  The p-group
#     version will have much weirder embedding maps.
        [IsPcGroup, IsPcGroup, IsGroupHomomorphism, IsListOrCollection],
        function ( N, Q, act, dom )
   local  F, G, col, ord, emb, embQ, embN, i, j, k, l, off;

   ord := Flat( Concatenation( RelativeOrders( Pcgs( Q ) ),
       ListWithIdenticalEntries( Size( dom ),
       RelativeOrders( Pcgs( N ) ) ) ) );
   # F := FreeGroup(Length(ord)); # boring names
   F := FreeGroup( Concatenation( List( Pcgs( Q ), String ),
       Concatenation( List( dom, x -> List( Pcgs( N ),
       g -> Concatenation( String( g ), "_", String( x )
       ) ) ) ) ) );

   col := SingleCollector( F, ord );

   emb := function ( w, off )
       return Product( [ 1 .. Length( w ) ], i ->
           GeneratorsOfGroup( F )[(i + off)] ^ w[i], One( F ) );
   end;
   embQ := w -> emb( w, 0 );
   embN := function ( w, i )
       return emb( w, Length( Pcgs( Q ) ) + Length( Pcgs( N ) ) * (i - 1) );
   end;

   for j  in [ 1 .. Length( Pcgs( Q ) ) ]  do
       SetPower( col, j, embQ( ExponentsOfRelativePower( Pcgs( Q ), j ) ) );
       for k  in [ 1 .. j - 1 ]  do
           SetConjugate( col, j, k, embQ( ExponentsOfConjugate( Pcgs( N ), j, k ) ) );
       od;
   od;
   
   for i  in [ 1 .. Length( dom ) ]  do
       off := Length( Pcgs( Q ) ) + Length( Pcgs( N ) ) * (i - 1);
       for j  in [ 1 .. Length( Pcgs( N ) ) ]  do
           SetPower( col, off + j, embN( ExponentsOfRelativePower( Pcgs( N ), j ), i ) );
           for k  in [ 1 .. j - 1 ]  do
               SetConjugate( col, off + j, off + k, embN( ExponentsOfConjugate( Pcgs( N ), j, k ), i ) );
           od;
       od;
   od;

   for i  in [ 1 .. Length( dom ) ]  do
       off := Length( Pcgs( Q ) ) + Length( Pcgs( N ) ) * (i - 1);
       for j  in [ 1 .. Length( Pcgs( Q ) ) ]  do
           k := dom[i]^Image(act,Pcgs( Q )[j]);
# WHY???           if i = k then continue; fi;
           for l  in [ 1 .. Length( Pcgs( N ) ) ]  do
               SetConjugate( col, off + l, j, F.(Length( Pcgs( Q ) ) + Length( Pcgs( N ) ) * (k - 1) + l) );
           od;
       od;
   od;

   G := GroupByRwsNC( col );

   # Setup embeddings
   SetWreathProductInfo( G, rec(l := Length(dom), q := Length(Pcgs(Q)),
       n := Length(Pcgs(N)), N := N, Q := Q, dom := dom, embeddings := []) );

   return G;
end);

InstallMethod(Embedding, "for pc wreath product",
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
            info.embeddings[i] := GroupHomomorphismByImagesNC(info.N,W,
                GeneratorsOfGroup(info.N), List(GeneratorsOfGroup(info.N),
    x->PcElementByExponentsNC(Pcgs(W), FilledIn(ExponentsOfPcElement(
        Pcgs(info.N),x),info.q+(i-1)*info.n,info.q+info.l*info.n))));
        elif i=info.l+1 then
            info.embeddings[i] := GroupHomomorphismByImagesNC(info.Q,W,
                GeneratorsOfGroup(info.Q), List(GeneratorsOfGroup(info.Q),
    x->PcElementByExponents(Pcgs(W), FilledIn(ExponentsOfPcElement(
        Pcgs(info.Q),x),0,info.q+info.l*info.n))));
        else
            return fail;
        fi;
        SetIsInjective(info.embeddings[i],true);
    fi;
    return info.embeddings[i];
end);
