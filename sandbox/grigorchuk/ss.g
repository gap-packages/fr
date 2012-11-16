#
# specialized code for the thinned group ring of the grigorchuk group
#

F := GF(2);
R := FreeAssociativeAlgebraWithOne(F,"A","B","C","D");
M := MatrixAlgebra(R,2);
R := LeftActingDomain(M);
A := R.1; B := R.2; C := R.3; D := R.4;

# phi := AlgebraWithOneHomomorphismByImages(R,M,[A,B,C,D],
#   [[[One(R),One(R)],[One(R),One(R)]],[[A,Zero(R)],[Zero(R),C]],
#    [[A,Zero(R)],[Zero(R),D]],[[Zero(R),Zero(R)],[Zero(R),B]]]);

################################################################
# the support of the expression x
#
Support := function(x)
  local e;
  e := CoefficientsAndMagmaElements(x);
  return List([1..Length(e)/2],n->e[2*n-1]);
end;

################################################################
# reduce using DAD=AA={B,C,D}{B,C,D}=0
#
WordCollapses := function(w)
  local e, i;
  e := ExtRepOfObj(w);
  for i in [1,3..Length(e)-1] do
    if e[i+1]>=2 then return true; fi;
    if i<=Length(e)-2 and e[i]>1 and e[i+2]>1 then return true; fi;
    if i<=Length(e)-4 and e[i]=4 and e[i+4]=4 then return true; fi;
  od;
  return false;
end;

mA := CoefficientsAndMagmaElements(A)[1];
mB := CoefficientsAndMagmaElements(B)[1];
mC := CoefficientsAndMagmaElements(C)[1];
mD := CoefficientsAndMagmaElements(D)[1];
mOne := CoefficientsAndMagmaElements(One(R))[1];

#
# convert a list (= external representation) to an element of same type as x
#
ListToMagma := function(l,x)
  local F;
  F := FamilyObj(x);
  return Objectify(F!.defaultType,
    NormalizedElementOfMagmaRingModuloRelations(F, [ ZeroCoefficient(x), l ]));
end;

Reduce := function(x)
  local e, i;
  if IsList(x) then return List(x,Reduce); fi;

  e := ShallowCopy(CoefficientsAndMagmaElements(x));
  for i in [1,3..Length(e)-1] do
    if WordCollapses(e[i]) then e[i+1] := ZeroCoefficient(x); fi;
  od;
  e := ListToMagma(FMRRemoveZero(e,ZeroCoefficient(x)),x);
  if Length(IntersectionSet(Support(e),[mB,mC,mD]))>=2 then
    e := e+B+C+D;
  fi;
  return e;
end;

################################################################
# the augmentation map A --> GF(2)
#
Epsilon := function(x)
  local e, i, z;
  if IsList(x) then return Epsilon(x[1][1])+Epsilon(x[1][2]); fi;

  if mOne in Support(x) then return Z(2)^0; else return 0*Z(2); fi;
end;

################################################################
# the map A --> M_2(A)
#
#b1000 := Basis(M)[1];
#b0100 := Basis(M)[2];
#b0010 := Basis(M)[3];
#b0001 := Basis(M)[4];
#mats := [b1000+b0100+b0010+b0001,b1000*A+b0001*C,b1000*A+b0001*D,b0001*B];

Word2Matrix := function(w)
#  local e, i, m;
#
#  e := ExtRepOfObj(w);
#  m := One(M);
#  for i in [1,3..Length(e)-1] do
#    if e[i+1]>=2 then return Zero(M); else m := m * mats[e[i]]; fi;
#  od;
#  return m;

  local e, i, v, pre, post, s, z;
  e := ExtRepOfObj(w);
  v := [[],[]];
  z := [false,false];
  pre := false; post := false;
  s := 1;
  for i in [1,3..Length(e)-1] do
    if e[i+1]>=2 then return Zero(M);
    elif e[i]=1 then
      if i=1 then pre := true; fi;
      s := 3-s;
      if i>1 and i=Length(e)-1 then post := true; fi;
    elif e[i]=2 then
      Append(v[s],[1,1]); Append(v[3-s],[3,1]);
    elif e[i]=3 then
      Append(v[s],[1,1]); Append(v[3-s],[4,1]);
    elif e[i]=4 then
      z[s] := true; Append(v[3-s],[2,1]);
    fi;
  od;
  if z[1] then v[1] := Zero(R); else v[1] := One(R)*ObjByExtRep(FamilyObj(w),v[1]); fi;
  if z[2] then v[2] := Zero(R); else v[2] := One(R)*ObjByExtRep(FamilyObj(w),v[2]); fi;
  if s=1 then e := [[v[1],Zero(R)],[Zero(R),v[2]]];
  else e := [[Zero(R),v[1]],[v[2],Zero(R)]];
  fi;
  if pre then e := Reduce([[One(R),One(R)],[One(R),One(R)]]*e); fi;
  if post then e := Reduce(e*[[One(R),One(R)],[One(R),One(R)]]); fi;
  return e;
end;
A2Matrix := function(x)
  local e, i, m;

  if IsList(x) then return List(x,A2Matrix); fi;
  e := CoefficientsAndMagmaElements(x);
  m := Zero(M);
  for i in [1,3..Length(e)-1] do
    m := m+e[i+1]*Word2Matrix(e[i]);
  od;
  return Reduce(m);
end;

################################################################
# returns true iff x is the zero element of A
#
FixIsZero := x-> (IsList(x) and ForAll(x,FixIsZero)) or (not IsList(x) and IsZero(x));
AIsZero := function(x)
  local z;
  if not IsZero(Epsilon(x)) then return false; fi;
  z := Reduce(x);
  return FixIsZero(z) or ForAll(Flat(A2Matrix(z)),AIsZero);
end;
AIsOne := x->AIsZero(A+One(R));

################################################################
# order of nilpotent element, rounded to upper power of 2
#
AOrder := function(x)
  local z;
  if not IsZero(Epsilon(x)) then return infinity;
  elif AIsZero(x) then return 1;
  else return 2*AOrder(Reduce(x^2));
  fi;
end;

################################################################
# basis of the augmentation ideal
#
BasisOfDelta := function(n)
  if n=0 then return [""];
  elif n=1 then return ["A","C","D"];
  else return Union(List(BasisOfDelta(n-1), function(x)
    local i;
    if x[1]='A' and 'C' in x{[2,6..4*QuoInt(Length(x)+2,4)-2]} then
      return [Concatenation("B",x)];
    elif x[1]='A' then return [Concatenation("B",x),Concatenation("C",x)];
    else return [Concatenation("A",x)];
    fi; end));
  fi;
end;

################################################################
# convert a word in ABCD to an element of R
#
W2R := function(w)
  return Product(List(ShallowCopy(w),function(x)
    if x='A' then return A;
    elif x='B' then return B;
    elif x='C' then return C;
    elif x='D' then return D;
    fi; end));
end;

################################################################
# convert a 2x2 matrix back to an element of A (if possible)
#
#sigma := AlgebraWithOneHomomorphismByImages(R,R,[A,B,C,D],[A*C*A,D,B,C]);
#
# * given x in A, M0001 computes z in A with A2Matrix(z) = [[*,0],[0,x]]
#   and * in Subring(A,[a,d])
# * given x,y in A, M1001 computes z in A with A2Matrix(z) = [[x,0],[0,y]],
#   or fail if there is no such z
# * given x in M2[A], Matrix2A computes z with A2Matrix(z) = x
#   or fail if there is no such z
#
Word2A := function(w)
  local e, f, i, pre, post;

  e := ExtRepOfObj(w);
  f := [];
  pre := One(R);
  post := One(R);
  for i in [1,3..Length(e)-1] do
    if e[i+1]>=2 then return Zero(R); fi;
    if e[i]=1 then
      if i=1 then pre := A-One(R); else Append(f,[1,1]); fi;
      Append(f,[3,1]);
      if i=Length(e)-1 then post := A-One(R); else Append(f,[1,1]); fi;
    elif e[i]=2 then Append(f,[4,1]);
    elif e[i]=3 then Append(f,[2,1]);
    elif e[i]=4 then Append(f,[3,1]);
    fi;
  od;
  f := pre*ObjByExtRep(FamilyObj(w),f)*post;
  if Length(e)=6 and e[1]>=3 and e[5]>=3 then
    return f+C*A*C*A*C;
  else return f;
  fi;
end;

M0001 := function(x)
  local e, i, z;

  e := CoefficientsAndMagmaElements(x);
  z := Zero(R);
  for i in [1,3..Length(e)-1] do z := z + e[i+1]*Word2A(e[i]); od;
  return z;
end;
M1001 := function(x,y)
  local xx, z;

  xx := Reduce(x);
  z := M0001(xx);
  z := Reduce((A-One(R))*z*(A-One(R))+M0001(Reduce(y-A2Matrix(z)[1][1])));
  if A2Matrix(z)[1][1]=xx then return z; else return fail; fi;
end;
Matrix2A := function(x)
  local i, j;

  if IsList(x[1][1]) then
    return [[Matrix2A(x[1][1]),Matrix2A(x[1][2])],
            [Matrix2A(x[2][1]),Matrix2A(x[2][2])]];
  else
    i := M1001(x[1][1],x[2][2]);
    j := M1001(x[1][2],x[2][1]);
    if i=fail or j=fail then return fail;
    else return Reduce(i+j*(A-One(R)));
    fi;
  fi;
end;

################################################################
# heavy duty reduction. build up a larger matrix, reduce, and re-build
# a ring element.
#
RR := function(x,level)
  if level=0 then return Reduce(x); else
    return Matrix2A(RR(A2Matrix(x),level-1));
  fi;
end;

ReduceALot := function(x)
  if IsList(x) then return List(x,ReduceALot);
  elif ForAll(Support(x),w->Length(w)<=2) then return Reduce(x);
  else return Matrix2A(ReduceALot(A2Matrix(x)));
  fi;
end;

################################################################
################################################################
################################################################

AttemptD := NewDictionary(A,false);
InverseD := NewDictionary(A,true);
IC := NewInfoClass("TraceInverse");
SetInfoLevel(IC,1);
#
# attempt to inverse an element in A
#
AInverse := function(x)
  local a, b, c, d, m, t, u, v;

  if IsZero(Epsilon(x)) then return fail; fi;
  t := ReduceALot(x);
  if IsOne(t) then return t; fi;
  if KnowsDictionary(AttemptD,t) then
    return LookupDictionary(InverseD,t);
  fi;
  AddDictionary(AttemptD,t);
  if Size(Support(t))=2 then
    u := One(R); v := Zero(R);
    Info(IC,2,"Power-invert element ",Support(t));
    while not AIsZero(u) do v := v+u; u := u*(t - One(R)); od;
    Info(IC,2,"Power-inverse is ",Support(v));
    v := ReduceALot(v);
    AddDictionary(InverseD,t,v);
    return v;
  fi;

  Info(IC,1,"Invert element ",Support(t));

  m := A2Matrix(t);
  a := m[1][1]; b := m[1][2]; c := m[2][1]; d := m[2][2];
  Info(IC,2,"Convert to matrix ",[Support(a),Support(b),Support(c),Support(d)]);

  if IsZero(Epsilon(a)) then
    v := AInverse(c);
    Info(IC,2,"Inverted c as ",Support(v));
    u := AInverse(b-a*v*d);
    Info(IC,2,"Computed c' as ",Support(u));
    v := [[-v*d*u,AInverse(c-d*AInverse(b)*a)],[u,-u*a*v]];
  else
    v := AInverse(d);
    Info(IC,2,"Inverted d as ",Support(v));
    u := AInverse(a-b*v*c);
    Info(IC,2,"Computed a' as ",Support(u));
    v := [[u,-u*b*v],[-v*c*u,AInverse(d-c*AInverse(a)*b)]];
  fi;
  Info(IC,2,"Computed M' as ",v);
  v := ReduceALot(Matrix2A(v));
  AddDictionary(InverseD,t,v);
  Info(IC,2,"Return inverse as ",Support(v));
  return v;
end;

Check := function(x)
  local z;
  z := AInverse(x);
  return [AIsZero(x*z-one),AIsZero(z*x-one)];
end;
