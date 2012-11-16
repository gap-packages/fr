f := FreeGroup("a","b","c");
ffff := FreeGroup("a1","b1","c1","a2","b2","c2","a3","b3","c3","a4","b4","c4");
fgens := List([1..4],i->GeneratorsOfGroup(ffff){[3*i-2..3*i]});
ffff := ffff / Flat([
                List(fgens,x->List(x,x->x^2)),
                List(Combinations([1..4],2),p->List(Cartesian(fgens{p}),LeftNormedComm)),
                List(fgens,x->(x[1]*x[3])^4),
                List(fgens,x->(x[1]*x[2])^8),
                List(fgens,x->(x[2]*x[3])^8),
                List(fgens,x->(x[3]^x[1]*x[3]^x[2])^2),
                List(fgens,x->Comm(Comm(x[2],x[3]),Comm(x[2],x[3])^x[1]))
                ]);
fgens := List([1..4],i->GeneratorsOfGroup(ffff){[3*i-2..3*i]});
n := Subgroup(ffff,[fgens[1][1]*fgens[2][3]*fgens[3][1]*fgens[4][3],
             fgens[1][3]*fgens[2][1]*fgens[3][3]*fgens[4][1],
             fgens[1][2]*fgens[4][2],
             fgens[2][2]*fgens[3][2]]);


#g := FRGroup("a=(1,2)(3,4)(5,6)","b=<a,c,a,c,,c>","c=<a,d,,d,a,d>","d=<,b,a,b,a,b>");
#g := FRGroup("a=(1,2)(3,4)","b=<a,c,a,c>","c=<b,,,b>");