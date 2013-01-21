RequirePackage("treegp");
#G := MakeFrGroup("a=[,,,,,](1,2)(3,4)(5,6)","b=[a,c,a,c,,c]()","c=[a,d,,d,a,d]()","d=[,b,a,b,a,b]()");
#AssignGeneratorVariables(G);

G := MakeFrGroup("A=[,,,](1,2)(3,4)","B=[,B,A,B]()","C=[A,C,,C]()","D=[A,D,A,D]()");
AssignGeneratorVariables(G);
pi := Decomposition(G);

MakeReadWriteGlobal("Delta");
Delta := FreeGroup("e1","e2","f1","f2","g2","h1");
AssignGeneratorVariables(Delta);
x := FreeGroup("e","f","g","h"); DeltaDelta := DirectProduct(x,x);
dec := GroupHomomorphismByImages(Delta,DeltaDelta,[e1,e2,f1,f2,g2,h1],[x.1^Embedding(DeltaDelta,1),x.1^Embedding(DeltaDelta,2),
               x.2^Embedding(DeltaDelta,1),x.2^Embedding(DeltaDelta,2),
               x.3^Embedding(DeltaDelta,1)*x.3^Embedding(DeltaDelta,2),
               x.4^Embedding(DeltaDelta,1)/x.4^Embedding(DeltaDelta,2)]);

Gamma := FreeGroup("a","b","c","d","t","u","v");
AssignGeneratorVariables(Gamma);
epi := GroupHomomorphismByImages(Gamma,G,[a,b,c,d,t,u,v],[A,B,C,D,One(G),One(G),One(G)]);
inj := GroupHomomorphismByImages(Delta,Gamma,[e1,e2,f1,f2,g2,h1],[(c*c^a)^2,(b*b^a)^2,((c*c^a)^2)^b,((b*b^a)^2)^c,d*d^a,c*d^a*b]);

test := x->IsOne(x^(inj*epi));

e1ij := function(i,j) return g2^i*h1^j*e1*h1^-j*g2^-i; end;
f1ij := function(i,j) return g2^i*h1^j*f1*h1^-j*g2^-i; end;
g1ij := function(i,j) return g2^i*h1^j*(h1^-1*e1*f1^-1)^j*g2^-i; end;
h2i  := function(i)   return h1^-1*(e1*f1^-1*g2)^i*g2^-i; end;

R1 := function(i,j) return Comm(e1ij(i,j),e2); end;
R2 := function(i,j) return Comm(e1ij(i,j),f2); end;
R3 := function(i,j) return Comm(f1ij(i,j),e2); end;
R4 := function(i,j) return Comm(f1ij(i,j),f2); end;
R5 := function(i,j) return Comm(g1ij(i,j),e2); end;
R6 := function(i,j) return Comm(g1ij(i,j),f2); end;
R7 := function(i,j,n) return f1ij(i,j)^-n*e1ij(i,j-1)^n*g1ij(i-1,j)^-1*f1ij(i-1,j)^-n*g1ij(i-1,j)*e1ij(i,j)^n; end;
R8 := function(i,j,n) return f1ij(i,j)^-n*g1ij(i,j)*e1ij(i+1,j)^(n-1)*g1ij(i,j+1)^-1*f1ij(i,j+1)^(1-n)*e1ij(i,j)^n; end;
R9 := function(i,n) return f2^-n*h2i(i)^-1*e2^n*h2i(i)*g2^-1*f2^-n*g2*e2^n; end;
R10 := function(i,n) return f2^-n*g2*e2^(n-1)*h2i(i-1)*g2^-1*f2^(1-n)*h2i(i)^-1*e2^n; end;

imax := 5;
nmax := 4;

if false then
Q := Delta / Concatenation(List(Cartesian([-imax..imax],[-imax..imax]),p->R1(p[1],p[2])),
             List(Cartesian([-imax..imax],[-imax..imax]),p->R2(p[1],p[2])),
             List(Cartesian([-imax..imax],[-imax..imax]),p->R3(p[1],p[2])),
             List(Cartesian([-imax..imax],[-imax..imax]),p->R4(p[1],p[2])),
             List(Cartesian([-imax..imax],[-imax..imax]),p->R5(p[1],p[2])),
             List(Cartesian([-imax..imax],[-imax..imax]),p->R6(p[1],p[2])),
             List(Cartesian([-imax..imax],[-imax..imax],[-nmax..nmax]),p->R7(p[1],p[2],p[3])),
             List(Cartesian([-imax..imax],[-imax..imax],[-nmax..nmax]),p->R8(p[1],p[2],p[3])),
             List(Cartesian([-imax..imax],[-imax..imax]),p->R9(p[1],p[2])),
             List(Cartesian([-imax..imax],[-imax..imax]),p->R10(p[1],p[2])));
else
    Q := Gamma / Concatenation([a^2,b^2,c^2,d^2,b*c*d,a^t/b^a,a^u/c^a,a^v/d^a],
         List(Cartesian([b,c,d],[t,u,v]),p->Comm(p[1],p[2])),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R1(p[1],p[2])^inj),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R2(p[1],p[2])^inj),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R3(p[1],p[2])^inj),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R4(p[1],p[2])^inj),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R5(p[1],p[2])^inj),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R6(p[1],p[2])^inj),
         List(Cartesian([-imax..imax],[-imax..imax],[-nmax..nmax]),p->R7(p[1],p[2],p[3])^inj),
         List(Cartesian([-imax..imax],[-imax..imax],[-nmax..nmax]),p->R8(p[1],p[2],p[3])^inj),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R9(p[1],p[2])^inj),
         List(Cartesian([-imax..imax],[-imax..imax]),p->R10(p[1],p[2])^inj));
fi;

F00 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^0,a^0,a^4,b]()");
F01 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^0,a^1,a^3,b]()");
F02 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^0,a^2,a^2,b]()");
F03 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^0,a^3,a^1,b]()");
F04 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^0,a^4,a^0,b]()");
F10 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^1,a^0,a^3,b]()");
F11 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^1,a^1,a^2,b]()");
F12 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^1,a^2,a^1,b]()");
F13 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^1,a^3,a^0,b]()");
F14 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^1,a^4,a^4,b]()");
F20 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^2,a^0,a^2,b]()");
F21 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^2,a^1,a^1,b]()");
F22 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^2,a^2,a^0,b]()");
F23 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^2,a^3,a^4,b]()");
F24 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^2,a^4,a^3,b]()");
F30 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^3,a^0,a^1,b]()");
F31 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^3,a^1,a^0,b]()");
F32 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^3,a^2,a^4,b]()");
F33 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^3,a^3,a^3,b]()");
F34 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^3,a^4,a^2,b]()");
F40 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^4,a^0,a^0,b]()");
F41 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^4,a^1,a^4,b]()");
F42 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^4,a^2,a^3,b]()");
F43 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^4,a^3,a^2,b]()");
F44 := MakeFrGroup("a=[,,,,](1,2,3,4,5)","b=[a,a^4,a^4,a^1,b]()");
