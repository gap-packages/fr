dirprod := 
    function( list )
    local len, D, F, G, pcgsG, gensF, s, h, i, j, t,
          info, first, coll, orders;

    # Check the arguments.
    if ForAny( list, G -> not IsPcGroup( G ) ) then
      TryNextMethod();
    fi;
    if ForAll( list, IsTrivial ) then
      return list[1];
    fi;
    len := Sum( List( list, x -> Length( Pcgs( x ) ) ) );
    F   := FreeGroup(IsSyllableWordsFamily, len );
    orders := [];
    for G in list do
        Append( orders, RelativeOrders(Pcgs(G)) );
    od;
    if ForAll(orders,x->x=orders[1]) then
        coll := CombinatorialCollector(F,orders);
    else
	coll := SingleCollector(F,orders);
    fi;
    gensF := GeneratorsOfGroup( F );

    s := 0;
    first := [1];
    for G in list do
        pcgsG := Pcgs( G );
        len   := Length(pcgsG);
        for i in [1..len] do
	    exp := ExponentsOfRelativePower( pcgsG, i );
	    t := One( F );
            for h in [1..len] do
                t := t * gensF[s+h]^exp[h];
            od;
            SetPower(coll, s+i, t);
	    for j in [i+1..len] do
	        exp := ExponentsOfPcElement( pcgsG, pcgsG[j]^pcgsG[i] );
            	t := One( F );
		for h in [1..len] do
                    t := t * gensF[s+h]^exp[h];
            	od;
            	SetConjugate(coll, s+j, s+i, t);
	    od;
	od;
        s := s+len;
        Add( first, s+1 );
    od;

    # create direct product
    D := GroupByRws(coll);

    # create info
    info := rec( groups := list,
                 first  := first,
                 embeddings := [],
                 projections := [] );
    SetDirectProductInfo( D, info );
    return D;
end;
