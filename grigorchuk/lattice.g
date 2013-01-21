Read("TreeGroups.g");
n := 7;
G := []; G[n] := MakeGGroup(n);
a := G[n].1; b := G[n].2; c := G[n].3; d := G[n].4; x := Comm(a,b);
NextIndex := function(g)
  local gG, s, t;
  gG := CommutatorSubgroup(g,G[n]);
  if Index(g,gG)=2 then return [gG];
  else
    s := First(GeneratorsOfGroup(g),x->not x in gG);
    t := First(GeneratorsOfGroup(g),x->not x in gG and not x*s in gG);
    return [ClosureGroup(gG,s),ClosureGroup(gG,t),ClosureGroup(gG,s*t)];
  fi;
end;

S := [];
Add(S,Set([[a,b],[a,c],[a,d],[b,c],[a*b,c],[a*c,d],[a*d,b]],p->ClosureGroup(DerivedSubgroup(G),p)));
Add(S,Set([a,b,c,d,a*b,a*c,a*d],p->ClosureGroup(DerivedSubgroup(G),p)));
Add(S,Set([[Comm(a,c),c*x],[c],[x,c^a*d],[x,Comm(a,c)],[b],[Comm(a,d),x*d],[d,x^2]],p->NormalClosure(G,Subgroup(G,p))));
Add(S,Set([[Comm(a,c)],[x],[Comm(a,d),x^2],[d],[Comm(a,d),x^2*d]],p->NormalClosure(G,Subgroup(G,p))));
Add(S,Set([[Comm(d,b^a),x^2],[Comm(a,d)*x^2],[Comm(a,d)]],p->NormalClosure(G,Subgroup(G,p))));
while Size(S)<35 do Add(S,Union(List(S[Size(S)],NextIndex))); od;
