if false then
    G := FreeGroup("a","b","c","d");
    AssignGeneratorVariables(G);
    G := G / Union([a^2,b^2,c^2,d^2,b*c*d,Comm(a,c)^2,Comm(a,d)^2],List([0,2..6],i->Comm(b,b^(a*(c*a)^i))));
    AssignGeneratorVariables(G);
    pi := GroupHomomorphismByImages(G,SymmetricGroup(6),[a,b,c,d],[(1,2),(3,4),(5,6),(3,4)(5,6)]);
    K := Kernel(pi);
    P := PresentationSubgroup(G,K);
    TzPrint(P);

    if false then
        "u=[b,a], v=[c,a], w=[d,a]";
        "K = < u,v,w,x | x=uw, uu, (Vux)^2, (xuV)^2, (VVuxx)^2, (xxuVV)^2, ... >";
    fi;

    D := FreeGroup("s","t","u");
    AssignGeneratorVariables(D);
    D := D / Union([s^2,t^2],List([-2..2],n->(s^-n*u*t^n)^2));
    AssignGeneratorVariables(D);

    DD := FreeGroup("s1","t1","u1","s2","t2","u2");
    AssignGeneratorVariables(DD);
    DD := DD / Union(List(Cartesian([s1,t1,u1],[s2,t2,u2]),p->Comm(p[1],p[2])),
                List([-2..2],n->(s1^-n*u1*t1^n)^2),
                List([-2..2],n->(s2^-n*u2*t2^n)^2));
    AssignGeneratorVariables(DD);
    R := Subgroup(DD,[t1*u2,u1*t2,s1*u1*s2*u2]);

    GG := FreeGroup("a1","b1","c1","d1","a2","b2","c2","d2");
    AssignGeneratorVariables(GG);
    GG := GG / Union([a1^2,b1^2,c1^2,d1^2,b1*c1*d1],
                 [a2^2,b2^2,c2^2,d2^2,b2*c2*d2],
                 List([0,2..8],i->Comm(b1,b1^(a1*(c1*a1)^i))),
                 List([0,2..8],i->Comm(b2,b2^(a2*(c2*a2)^i))),
                 List(Cartesian([a1,b1,c1,d1],[a2,b2,c2,d2]),p->Comm(p[1],p[2])));
    AssignGeneratorVariables(GG);
    R := Subgroup(GG,[a1*a2,b1*c2,c1*b2,d1*d2]);

    pi := PqEpimorphism(GG : Prime := 2, ClassBound := 12);
    Q := Image(pi);
    R := Image(pi,R);
    R0 := Core(Q,R);
    S := R/R0;
    DerivedSeries(S);
fi;

RequirePackage("treegp");
G := MakeFrGroup("a=[,,,](1,2)(3,4)","b=[,b,a,b]()","c=[a,c,,c]()","d=[a,d,a,d]()");
AssignGeneratorVariables(G);
pi := Decomposition(G);

e2 := (c*c^a)^2;
