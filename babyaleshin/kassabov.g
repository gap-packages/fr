LoadPackage("fr");
g := PermGroup(BabyAleshinGroup,10);
x := Comm(g.1,g.2);
y := Comm(g.2,g.3);
z := Comm(g.3,g.1);
lcs := PCentralSeries(Group(x,y,z));

u := y*z;
v := y;
h := x*y*z;

extend := function(n,l)
    local r, i, j, k, gu, gv, gh;
    r := [];
    for l in l do
        for i in Cartesian([0,2^n],[0,2^n],[0,2^n]) do
            gu := x^(l[1][1]+i[1])*y^(l[1][2]+i[2])*z^(l[1][3]+i[3]);
            for j in Cartesian([0,2^n],[0,2^n],[0,2^n]) do
                gv := x^(l[2][1]+j[1])*y^(l[2][2]+j[2])*z^(l[2][3]+j[3]);
                for k in Cartesian([0,2^n],[0,2^n],[0,2^n]) do
                    gh := x^(l[3][1]+k[1])*y^(l[3][2]+k[2])*z^(l[3][3]+k[3]);
                    if Comm(gu,gv)/gh^4 in lcs[n+4] and
                       Comm(gu,gh)/gu^8 in lcs[n+4] and
                       Comm(gv,gh)*gv^8 in lcs[n+4] then
                        Add(r,[l[1]+i,l[2]+j,l[3]+k]);
                    fi;
                od;
            od;
        od;
    od;
    return r;
end;

next := function(g,l,n)
    return First(Cartesian(l[1]+[0,2^n],l[2]+[0,2^n],l[3]+[0,2^n]),
                   p->g/(u^p[1]*v^p[2]*h^p[3]) in lcs[n+2]);
end;

l0 := [[[1,0,0],[0,1,0],[1,1,1]]];

n := 9;
o := One(Integers mod 2^n);
G := Group([[1,8],[0,1]]*o,[[1,0],[1,1]]*o,[[5,0],[0,1/5]]*o);
pi := ActionHomomorphism(G,[1..2^n],function(p,x) local z; z := ZmodnZObj(p-1,2^n); z := (x[1][1]*z+x[2][1])/(x[1][2]*z+x[2][2]); return z![1]+1; end);
Gp := Image(pi);

n := 10;
g := PermGroup(BabyAleshinGroup,n);
x := g.1*g.2;
y := g.2*g.3;
z := g.3*g.1;
g := Group(x,y,z);

conj := RepresentativeAction(SymmetricGroup(2^n),y,Gp.2);
#IsomorphismGroups(g,G);
l := Gp^(conj^-1);
l := Filtered(Centralizer(SymmetricGroup(2^n),y),i->x^i in l and z^i in l);
