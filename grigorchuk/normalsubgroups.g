Pi := function(n)
  local x, t, pi, sigma;
  pi := ['B','a'];
  x := pi;
  sigma := function(x)
    if x='a' then return ['a','D','a'];
    elif x='B' then return 'D';
    elif x='C' then return 'B';
    elif x='D' then return 'C';
    fi;
  end;
  for t in [1..n] do
    x := Flat(List(x,sigma));
    Append(pi,x);
  od;
  return pi;
end;

Pi := Pi(5);

TwoFloor := x->2^LogInt(x,2);

type1 := function(C,A)
  return List(Filtered(Combinations([Int(TwoFloor(A)/2)..C-1]),
    function(B)
      return C<=TwoFloor(A) and A<=5*TwoFloor(C) and
	ForAll(B,x->TwoFloor(C)/2<=x and A<5*TwoFloor(x) and
		(TwoFloor(x)=TwoFloor(C) or TwoFloor(x)=TwoFloor(A)/2)) and
	ForAll(B,x->Pi{[x..C-2]}=Pi{[A..C-2+A-x]});
      end),B->[A+C+2,"W(k)*PROD W(l)^2,W(j)^2",A,B,C]);
end;

type2 := function(C,A)
  return List(Filtered(Combinations([1..C-1]),
    function(B)
      return C>TwoFloor(A) and C<=2*TwoFloor(A) and
	A>=3*TwoFloor(A)/2 and
	Size(B)>=1 and B[1]=A-TwoFloor(A) and
	ForAll(B{[2..Size(B)]},x->Pi{[x..C-2]}=Pi{[A..C-2+A-x]});
      end),B->[A+C+2,"W(k)*PROD W(l)^2,W(j)^2",A,B,C]);
end;

type3 := function(C)
  return List(Filtered(Combinations([1..C-1]),
    function(B)
      return C<=TwoFloor(C)*3/2 and
	Size(B)>=1 and B[Size(B)]-B[1]=C-TwoFloor(C) and
	TwoFloor(B[Size(B)])=TwoFloor(C) and
	ForAll(B{[1..Size(B)-1]},x->TwoFloor(x)=TwoFloor(C)/2) and
	ForAll(B{[1..Size(B)-1]},x->Pi{[x..TwoFloor(C)-2]}=Pi{[B[Size(B)]..B[Size(B)]-x+TwoFloor(C)-2]});
      end),B->[B[1]+C+2+TwoFloor(C),"W(infinity)*PROD W(l)^2,W(j)^2",infinity,B,C]);
end;

DescToGroup := function(G,l)
  local N;
  x := Comm(G.1,G.2);
  N := NormalClosure(G,Subgroup(G,[W(l[3])*Product(l[4],WW,()),WW(l[5])]));
  if Index(G,N)<>2^l[1] then Error("Bad Index in <desc>",l); fi;
  return N;
end;

if false then
  LogTo("x");
  Print(Union(List(Cartesian([1..8],[1..16]),p->type1(p[1],p[2]))));
  Print(Union(List(Cartesian([1..8],[1..16]),p->type2(p[1],p[2]))));
  Print(Union(List([1..8],p->type3(p))),"\n");
  LogTo();
fi;

#LogTo("y");
#l := List(Combinations([4..10]),x->NS(W(16)*Product(x,WW,())));
#List(l,A);
#LogTo();
