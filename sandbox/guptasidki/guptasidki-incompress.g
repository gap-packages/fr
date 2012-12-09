tst := function(S,n)
  local i, j, k, x, y;
#  S := [[t,t^2]];
  for i in [Length(S)+1..n] do
    Add(S,[]);
    for j in [[1,1],[1,-1],[-1,1],[-1,-1]] do for k in S[i-1] do
      x := t^j[1]*a^j[2]*k;
      y := List(Decomposition(x)[1],gnum);
      if Sum(y)=i and Product(y)=0 then
        Add(S[i],x);
      fi;
    od; od;
  od;
  return S;
end;
g0 := [a,a^-1,a^0];
g1 := [t,t^-1];
gnum := function(x)
  if x in g0 then return 0;
  elif x in g1 then return 1;
  else return Sum(List(Decomposition(x)[1],gnum)); fi;
end;
  
more := function(m,n)
  local i, j, k, r, x;
  r := 0;
  for i in s[m] do for j in [1,-1] do for k in s[n] do
    x := i*a^j*k;
    if gnum(x)=m+n then
      r := r+1;
    fi;
  od; od; od;
  return r;
end;
