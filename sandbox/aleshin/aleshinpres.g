exps := function(g,v,l)
  local e,i;
  e := List(v,x->0);
  for i in [1..Length(l)-1] do
    e := e+First(2^(i-1)*Cartesian(List(v,x->[0..1])),p->Product([1..Length(v)],j->v[j]^(e[j]+p[j]))/g in l[i+1]);
  od;
  return e;
end;

pi := Decomposition(G);

vanishes := function(g,n)
  local x;
  if n=0 then return true;
  else x := pi(g);
    return x[2]=() and ForAll(x[1],g->vanishes(g,n-1));
  fi;
end;
  
F := FreeGroup("e","f","g");
AssignGeneratorVariables(F);
n := 13;
alpha := 357; beta := 581; gamma := 60; delta := 141; epsilon := 1524;
F := F/ [e^(2^(n-1)),f^(2^(n-2)),g^(2^(n-1)),e^(2^(n-4))*g^(2^(n-4)),
  f^(2^(n-3))*g^(2^(n-2)),
  f*e / (e^alpha*f^beta*g^gamma),
  g*f / (e^gamma*f^beta*g^alpha),
  g*e / (e^delta*f^epsilon*g^delta)];
