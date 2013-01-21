MakeGPres := function(n)
    local G, a, b, c, d, sigma;
    G := FreeGroup("a","b","c","d");
    a := G.1; b := G.2; c := G.3; d := G.4;
    sigma := GroupHomomorphismByImages(G,G,[a,b,c,d],[c^a,d,b,c]);
    return G / Union([a^2,b^2,c^2,d^2,b*c*d],
                   List([0..n+1], i->Image(sigma^i,(a*d)^4)),
                   List([0..n], i->Image(sigma^i,(a*d*a*c*a*c)^4)));
end;

MakeGPresF := function(n)
    local G, a, b, c, d, sigma;
    G := FreeGroup("a","b","c","d");
    a := G.1; b := G.2; c := G.3; d := G.4;
    sigma := GroupHomomorphismByImages(G,G,[a,b,c,d],[c^a,d,b,c]);
    return G / Union([a^2,b^2,c^2,d^2,b*c*d],
                   List([0..n+1], i->Image(sigma^i,(a*d)^4)),
                   List([0..n], i->Image(sigma^i,(a*d*a*c*a*c)^4)),
		   [Image(sigma^(n-2),(a*d)^2),
		    Image(sigma^(n-1),(a*b)^2),
		    Image(sigma^(n-2),(a*c)^4)]);
end;

MakeSPres := function(n)
    local G, a, b, c, d, sigma;
    G := FreeGroup("a","b","c","d");
    a := G.1; b := G.2; c := G.3; d := G.4;
    sigma := GroupHomomorphismByImages(G,G,[a,b,c,d],[b^a,d,b,c]);
    return G / Union([a^2,b^2,c^2,d^2,Comm(b,c),Comm(b,d),Comm(c,d)],
                   List([0..n+1], i->Image(sigma^i,(a*c)^4)),
                   List([0..n+1], i->Image(sigma^i,(a*d)^4)),
                   List([0..n+1], i->Image(sigma^i,(a*c*a*d)^2)),
                   List([0..n+1], i->Image(sigma^i,(a*b)^8)),
                   List([0..n], i->Image(sigma^i,(a*b*a*b*a*c)^4)),
                   List([0..n], i->Image(sigma^i,(a*b*a*b*a*d)^4)),
                   List([0..n], i->Image(sigma^i,(a*b*a*b*a*c*a*b*a*b*a*d)^2)));
end;
