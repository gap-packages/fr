# Read("Grigorchuk.g");

p := 2;
F := GF(p);
G := FreeGroup("a","b","c","d");
ga := G.1; gb := G.2; gc := G.3; gd := G.4;
G := G / [ga^2,gb^2,gc^2,gd^2,gb*gc*gd,(ga*gd)^4];
ga := G.1; gb := G.2; gc := G.3; gd := G.4;

H := Subgroup(G,[gb,gc,gd,gb^ga,gc^ga,gd^ga]);
#phi := [GroupHomomorphismByImages(H,G,[ga^2,gb,gc,gd,gb^ga,gc^ga,gd^ga],
#  [Identity(G),ga,ga,Identity(G),gc,gd,gb]),
#        GroupHomomorphismByImages(H,G,[ga^2,gb,gc,gd,gb^ga,gc^ga,gd^ga],
#  [Identity(G),gc,gd,gb,ga,ga,Identity(G)])];

_phi := function(x,s0)
  local e, i, v, s, P;
  v := Identity(G);
  e := ExtRepOfObj(x);
  s := s0;
  P := [[,gc,gd,gb],[,ga,ga,Identity(G)]];
  for i in [1,3..Length(e)-1] do
    if IsOddInt(e[i+1]) then
      if e[i] = 1 then s := not s;
      elif s then v := v*P[1][e[i]];
      else v := v*P[2][e[i]];
      fi;
    fi;
  od;
  return v;
end;

#
# given a word w stabilizing the first level, the action on the two subtrees
# are Image(phi[i],w) for i=1,2
#
phi := [GroupHomomorphismByFunction(H,G,x->_phi(x,false)),
        GroupHomomorphismByFunction(H,G,x->_phi(x,true))];

A := GroupRing(F,G);
one := One(A); a := A.1; b := A.2; c := A.3; d := A.4;

#
# reduce a word in G following the relations g^2=1, bc=d and adada=dad
#
ReduceWord := function(x)
  local i, e, changed;

  e := ExtRepOfObj(x);
  repeat
    changed := false;
				# reduce s^n, with n outside [1..p-1]
    for i in [1,3..Length(e)-1] do
      if e[i+1] >= p or e[i+1] <= 0 then
        changed := true;
        e[i+1] := e[i+1] mod p;
	if e[i+1]=0 then RemoveElmList(e,i); RemoveElmList(e,i); break; fi;
      fi;
    od;
				# reduce s^ns^m
    for i in [1,3..Length(e)-3] do
      if e[i] = e[i+2] then
        changed := true;
        e[i+1] := e[i+1] + e[i+3];
        RemoveElmList(e,i+2); RemoveElmList(e,i+2);
        break;
      fi;
    od;
				# reduce b^nc^n=d^n etc
    for i in [1,3..Length(e)-3] do
      if e[i] >= 2 and e[i+2] >= 2 and e[i] <> e[i+2] and e[i+1] = e[i+3] then
        changed := true;
        e[i] := 9 - e[i] - e[i+2];
        RemoveElmList(e,i+2); RemoveElmList(e,i+2);
        break;
      fi;
    od;
				# reduce dada -> adad
    for i in [1,3..Length(e)-7] do
      if e{[i..i+7]} = [4,1,1,1,4,1,1,1] then
        changed := true;
        e[i] := 1; e[i+2] := 4; e[i+4] := 1; e[i+6] := 4;
        break;
      fi;
    od;
				# reduce {b,c}adad$ -> {b,c}dada$
    for i in [Maximum(1,Length(e)-9)..Length(e)-9] do
      if e[i] in [2,3] and e{[i+1..i+9]} = [1,1,1,4,1,1,1,4,1] then
        changed := true;
        e[i+2] := 4; e[i+4] := 1; e[i+6] := 4; e[i+8] := 1;
        break;
      fi;
    od;
  until not changed;
  return AssocWord(TypeObj(x)![AWP_PURE_TYPE], e{[1..Length(e)]});
end;

#
# convert a list (= external representation) to an element of same type as x
#
ListToMagma := function(l,x)
  local F;
  F := FamilyObj(x);
  return Objectify(F!.defaultType,
    NormalizedElementOfMagmaRingModuloRelations(F, [ ZeroCoefficient(x), l ]));
end;

#
# the support of the expression x
#
Support := function(x)
  local e;
  e := CoefficientsAndMagmaElements(x);
  return List([1..Length(e)/2],n->e[2*n-1]);
end;

#
# reduce an element of A by applying reductions in G and collapsing
# coefficients
#
Reduce := function(x)
  local i, e, v, w;
  if IsList(x) then return List(x,Reduce); fi;
  if x in G then return ReduceWord(x); fi;

  e := CoefficientsAndMagmaElements(x);
  v := List([1..Length(e)/2],n->e[2*n-1]);
  w := List([1..Length(e)/2],n->e[2*n]);
  for i in [1..Length(v)] do
    v[i] := ReduceWord(v[i]);
  od;
  SortParallel(v,w);
  for i in [1..Length(v)-1] do
    if v[i] = v[i+1] then
      w[i+1] := w[i+1] + w[i];
      w[i] := ZeroCoefficient( x );
    fi;
  od;
  e := ListToMagma(FMRRemoveZero(List([1..2*Length(v)],
    function(n) if IsOddInt(n) then return v[(n+1)/2];
      else return w[n/2]; fi; end),
    ZeroCoefficient(x)),x);
  if Length(IntersectionSet(Support(e),[gb,gc,gd]))>=2 then
    e := e+one+b+c+d;
  fi;
  if Length(IntersectionSet(Support(e),[ga*gd,gd*ga,gd*ga*gd]))>=2 then
    e := e+a+a*d+d*a+d*a*d;
  fi;
  if Length(IntersectionSet(Support(e),[gd,ga*gd*ga,ga*gd*ga*gd]))>=2 then
    e := e+one+d+a*d*a+a*d*a*d;
  fi;
  return e;
end;

#
# the augmentation map A --> GF(2)
#
Epsilon := function(x)
  local e, i, z;

  e := CoefficientsAndMagmaElements(x);
  z := ZeroCoefficient(x);
  for i in [2,4..Length(e)] do z := z+e[i]; od;
  return z;
end;

#
# returns true iff x is the zero element of A
#
AIsZero := function(x)
  local z;
  z := Reduce(x);
  if IsZero(z) then return true; fi;
  z := Flat(A2Matrix(z));
  if ForAll(z,x->IsZero(Epsilon(x))) then return ForAll(z,AIsZero);
  else return false;
  fi;
end;

#
# convert an element of A to a 2x2 matrix over A
#
A2Matrix := function(x)
  local e, i, a, b, c, d;
  if IsList(x) then return List(x,A2Matrix); fi;

  e := CoefficientsAndMagmaElements(x);
  a := []; b := []; c := []; d := [];

  for i in [1,3..Length(e)-1] do
    if e[i] in H then
      Append(a,[Image(phi[1],e[i]),e[i+1]]);
      Append(d,[Image(phi[2],e[i]),e[i+1]]);
    else
      Append(b,[Image(phi[1],e[i]*ga),e[i+1]]);
      Append(c,[Image(phi[2],e[i]*ga),e[i+1]]);
    fi;
  od;
  return Reduce([[ListToMagma(a,x),ListToMagma(b,x)],
                 [ListToMagma(c,x),ListToMagma(d,x)]]);
end;

#
# given words x,y, s(x,y) is a word stabilizing the first level with
# Image(phi[i],s(x,y)) = x,y
# this really should be a group homomorphism
#
sigma := GroupHomomorphismByImages(G,H,[ga,gb,gc,gd],[gc^ga,gd,gb,gc]);
s := function(x,y)
  local s;
  s := Image(sigma,x);
  return Reduce(s^ga*Image(sigma,Image(phi[1],s^-1)*y));
end;

#
# @ given x in A, sigmaA computes z in A with A2Matrix(z) = [[*,0],[0,x]]
#   and * in Subring(A,[a,d])
# @ given x,y in A, sA computes z in A with A2Matrix(z) = [[x,0],[0,y]]
# @ given x in M2[A], Matrix2A computes z with A2Matrix(z) = x
#
sigmaA := function(x)
  local e, i, z;

  e := CoefficientsAndMagmaElements(x);
  z := [];
  for i in [1,3..Length(e)-1] do Append(z,[Image(sigma,e[i]),e[i+1]]); od;
  return ListToMagma(z,x);
end;
sA := function(x,y)
  local z;
  z := sigmaA(x);
  z := a*z*a+sigmaA(y-A2Matrix(z)[1][1]);
# a good first try... check to see if we added some (1+a)(1+a^d)
  if not AIsZero(A2Matrix(z)[1][1]-x) then
    z := Reduce(z-a*sigmaA(Reduce(A2Matrix(z)[1][1]-x))*a);
  fi;
  if not AIsZero(A2Matrix(z)[2][2]-y) then # something f@cked up
    Error("Hello, my name's Goofy!");
  fi;
  return z;
end;
Matrix2A := function(x)
  if IsList(x[1][1]) then
    return [[Matrix2A(x[1][1]),Matrix2A(x[1][2])],
            [Matrix2A(x[2][1]),Matrix2A(x[2][2])]];
  else
    return Reduce(sA(x[1][1],x[2][2])+sA(x[1][2],x[2][1])*a);
  fi;
end;

#
# heavy duty reduction. build up a larger matrix, reduce, and re-build
# a ring element.
#
RR := function(x,level)
  if level=0 then return Reduce(x); else
    return Matrix2A(RR(A2Matrix(x),level-1));
  fi;
end;

ID := NewDictionary(false,true,A);
IC := NewInfoClass("TraceInverse");
#
# attempt to inverse an element in A
#
AInverse := function(x)
  local a, b, c, d, t, u;

  if IsZero(Epsilon(x)) then return fail; fi;
  t := Reduce(x);
  if Length(Support(t))=1 then return Reduce(Inverse(t)); fi;

  Info(IC,1,"Invert element with support ",Length(Support(t)));
### UGLY! we're probably stuck in some sort of recursion
  if Length(Support(t))>99 then Error("I'm probably stuck!"); fi;


  t := A2Matrix(t);
  a := t[1][1]; b := t[1][2]; c := t[2][1]; d := t[2][2];

  if IsZero(Epsilon(a)) then
    t := AInverse(c);
    u := AInverse(b-a*t*d);
    t := [[-t*d*u,AInverse(c-d*AInverse(b)*a)],[u,-u*a*t]];
  else
    t := AInverse(d);
    u := AInverse(a-b*t*c);
    t := [[u,-u*b*t],[-t*c*u,AInverse(d-c*AInverse(a)*b)]];
  fi;
AddDictionary(ID,x,t);
  return Matrix2A(t);
end;

Check := function(x)
  local z;
  z := AInverse(x);
  return [AIsZero(x*z-one),AIsZero(z*x-one)];
end;
