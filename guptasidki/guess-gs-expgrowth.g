AssignGeneratorVariables(GuptaSidkiGroup);
U := [[a^-1]];
V := [[a]];
for n in [2..5] do
    Add(U,Concatenation(List(Cartesian(U[n-1],U[n-1]),p->p[1]^-1*t^-1*p[2]),
            List(Cartesian(V[n-1],V[n-1]),p->p[1]^-1*t*p[2])));
    Add(V,Concatenation(List(Cartesian(U[n-1],U[n-1]),p->p[1]^-1*t*p[2]),
            List(Cartesian(V[n-1],V[n-1]),p->p[1]^-1*t^-1*p[2])));
od;

    