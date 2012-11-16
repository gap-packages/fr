F := GF(2);
R := FreeAssociativeAlgebraWithOne(F,"A","B","C","D");
A := R.1; B := R.2; C := R.3; D := R.4;
zero := Zero(R); one := One(R);

################################################################
# reduce an algebra element using the relations
# AA = BB = CC = DD = BC = CB = BD = DB = CD = DC = B+C+D = 0
# D*****D = 0 (4n+1 stars)
################################################################
IsTrivialGGroupWord := function(x)
  local e, i, Dpos;

  e := ExtRepOfObj(x);
  Dpos := -1;

  for i in [2,4..Length(e)] do
    if e[i]>=2 then return true; fi;
    if i>=4 and e[i-1]>=2 and e[i-3]>=2 then return true; fi;
    if e[i-1]=4 then
      if Dpos>=0 and RemInt(i-Dpos,8)=4 then return true; fi;
      Dpos := i;
    fi;
  od;
  return false;
end;

# express x as a sum of words in A,B,D only, with no identical consecutive
# letters and no D*****D (4n+1 stars).
# it returning the set of these words.
# (this extends IsTrivialGGroupWord).
_CGGW := function(w,Dpos,x,result)
  local Cpos, Cmod;

  Cpos := First([1,3..Length(w)-1],i->w[i]=3);
  if Cpos = fail then
    Add(result,ObjByExtRep(FamilyObj(x),w));
  else
    Cmod := RemInt(Cpos,8);
    if Size(Union(Dpos,[Cmod]))>1 then
      w[Cpos] := 2;
      _CGGW(w,Dpos,x,result);
      w[Cpos] := 3;
    else
      w[Cpos] := 2;
      _CGGW(w,Dpos,x,result);
      w[Cpos] := 4;
      _CGGW(w,Union(Dpos,[Cmod]),x,result);
      w[Cpos] := 3;
    fi;
  fi;
end;

CanonizeGGroupWord := function(x)
  local e, i, Dpos, result;

  e := ExtRepOfObj(x);
  if ForAny([2,4..Length(e)],i->e[i]>=2) then return []; fi;
  if Size(Set(List(Filtered([1,3..Length(e)-1],i->e[i]<>1),i->RemInt(i,4))))>1 then return []; fi;
  Dpos := Set(List(Filtered([1,3..Length(e)-1],i->e[i]=4),i->RemInt(i,8)));
  if Size(Dpos)>1 then return []; fi;

  result := [];
  _CGGW(ShallowCopy(e),Dpos,x,result);
  return result;
end;

BCD := List([B,C,D],x->CoefficientsAndMagmaElements(x)[1]);

if false then
ThinnedReduce := function(x)
  local e, elm, coeff, i, BCDpos;

  if IsList(x) then return List(x,ThinnedReduce); fi;

  e := CoefficientsAndMagmaElements(x);
  elm := e{[1,3..Length(e)-1]};
  coeff := e{[2,4..Length(e)]};
  BCDpos := -1;

  for i in [1..Length(elm)] do
    if IsTrivialGGroupWord(elm[i]) then coeff[i] := ZeroCoefficient(x); fi;
    if elm[i] in BCD then
      if BCDpos=-1 then BCDpos := i;
      else coeff[BCDpos] := ZeroCoefficient(x);
        elm[i] := BCD[6-Position(BCD,elm[i])-Position(BCD,elm[BCDpos])];
        BCDpos := -1;
      fi;
    fi;
  od;

  return ElementOfMagmaRing(FamilyObj(x),ZeroCoefficient(x),coeff,elm);
end;
else
ThinnedReduce := function(x)
  local e;
  if IsList(x) then return List(x,ThinnedReduce); fi;

  e := CoefficientsAndMagmaElements(x);
  e := Flat(List(e{[1,3..Length(e)-1]},CanonizeGGroupWord));

  return ElementOfMagmaRing(FamilyObj(x),ZeroCoefficient(x),List(e,x->Z(2)),e);
end;
fi;
  
################################################################
# pretty-print
################################################################
GGroupWord2String := function(x)
  local e, i, s, t;

  e := ExtRepOfObj(x);
  s := "";
  t := ["A","B","C","D"];
  for i in [1,3..Length(e)-1] do
    s := Concatenation(s,t[e[i]]);
  od;

  if Length(e)=0 then return "1"; else return s; fi;
end;

Thinned2String := function(x)
  local e, i, s;

  e := CoefficientsAndMagmaElements(ThinnedReduce(x));
  s := "";

  for i in [1,3..Length(e)-1] do
    if i<>1 then s := Concatenation(s,"+"); fi;
    s := Concatenation(s,GGroupWord2String(e[i]));
  od;

  return s;
end;

################################################################
# the augmentation map A --> GF(2)
#
# the support of the expression x
################################################################
Epsilon := function(x)
  local e, i;

  if IsList(x) then return Sum(List(x,Epsilon)); fi;

  e := CoefficientsAndMagmaElements(x);
  for i in [2,4..Length(e)] do if IsOne(e[i-1]) then return e[i]; fi; od;
  return Zero(F);
end;

ThinnedSupport := function(x)
  local e;
  e := CoefficientsAndMagmaElements(x);
  return List([1..Length(e)/2],n->e[2*n-1]);
end;

################################################################
# convert an element of R to a 2x2 matrix over R
################################################################
#psi := AlgebraHomomorphismByImages(R,MatrixAlgebra(R,2),[A,B,C,D],
#  [[[one,one],[one,one]],
#   [[A,zero],[zero,C]],
#   [[A,zero],[zero,D]],
#   [[zero,zero],[zero,B]]]);

GGroupWord2Matrix := function(x)
  local e, i, u, pos;

  e := ExtRepOfObj(x);
  u := [one,one];
  pos := 1;

  for i in [1,3..Length(e)-1] do
    if e[i]=1 then pos := 3-pos;
    elif e[i]=2 then u[pos] := u[pos]*A; u[3-pos] := u[3-pos]*C;
    elif e[i]=3 then u[pos] := u[pos]*A; u[3-pos] := u[3-pos]*D;
    elif e[i]=4 then u[pos] := zero; u[3-pos] := u[3-pos]*B;
    fi;
  od;
  if pos=1 then u := [[u[1],zero],[zero,u[2]]];
  else u := [[zero,u[1]],[u[2],zero]];
  fi;
  if Length(e)>=2 and e[1]=1 then
     u := [[one,one],[one,one]]*u; fi;
  if Length(e)>=4 and e[Length(e)-1]=1 then
     u := u*[[one,one],[one,one]]; fi;

  return u;
end;

Thinned2Matrix := function(x)
  local e, i, m;
  if IsList(x) then return List(x,Thinned2Matrix); fi;

  e := CoefficientsAndMagmaElements(ThinnedReduce(x));
  m := [[zero,zero],[zero,zero]];

  for i in [1,3..Length(e)-1] do
    m := m + e[i+1]*GGroupWord2Matrix(e[i]);
  od;
  return ThinnedReduce(m);
end;

################################################################
# returns true iff x is the zero element of A
################################################################
FixIsZero := x-> (IsList(x) and ForAll(x,FixIsZero)) or (not IsList(x) and IsZero(x));

ThinnedIsZero := function(x)
  local z;
  if not IsZero(Epsilon(x)) then return false; fi;
  z := ThinnedReduce(x);
  return FixIsZero(z) or ForAll(Flat(Thinned2Matrix(z)),ThinnedIsZero);
end;

################################################################
# convert a 2x2 matrix over R to a matrix over R.
# this is the inverse of Thinned2Matrix
################################################################
GGroupSigma := function(x)
  local e, i, y;

  e := ExtRepOfObj(x);
  y := one;

  for i in [1,3..Length(e)-1] do
    if e[i]=1 then
      if i=1 then y := y*(A+one); else y := y*A; fi;
      y := y*C;
      if i=Length(e)-1 then y := y*(A+one); else y := y*A; fi;
    elif e[i]=2 then y := y*D;
    elif e[i]=3 then y := y*B;
    elif e[i]=4 then y := y*C;
    fi;
  od;
  return y;
end;

# this returns y=[[*,0],[0,x]] with * in <A,D>
_M2T := function(x)
  local e, i, y;

  e := CoefficientsAndMagmaElements(x);
  y := zero;

  for i in [2,4..Length(e)] do
    y := y+e[i]*GGroupSigma(e[i-1]);
  od;
  return y;
end;

# this returns z=[[x,0],[0,y]]
__M2T := function(x,y)
  local z, mz;

  z := _M2T(x);
  z := (A+one)*z*(A+one) + _M2T(y-Thinned2Matrix(z)[1][1]);
  mz := Thinned2Matrix(z);
  if IsZero(mz[1][1]-x) then
  elif IsZero(mz[1][1]-x-A*D*A) then
    z := z + C*A*C*A*C;
  else
    Error("These elements do not lift to R: ", x,y);
  fi;
  return z;
end;

Matrix2Thinned := function(x)
  local rx, z;

  if IsList(x[1][1]) then
    return [[Matrix2Thinned(x[1][1]),Matrix2Thinned(x[1][2])],
            [Matrix2Thinned(x[2][1]),Matrix2Thinned(x[2][2])]];
  else
    rx := ThinnedReduce(x);
    z := ThinnedReduce(__M2T(rx[1][1],rx[2][2])+__M2T(rx[1][2],rx[2][1])*(A+one));
    return z;
  fi;
end;

# heavy duty reduction. build up a larger matrix, reduce, and re-build
# a ring element.
RecursiveReduce := function(x,level)
  if level=0 then return ThinnedReduce(x);
  else
    return Matrix2Thinned(RecursiveReduce(Thinned2Matrix(x),level-1));
  fi;
end;

SuperReduce := function(x)
  if x=zero then return x; fi;
  return RecursiveReduce(x,LogInt(Maximum(List(ThinnedSupport(x),Length)),2)+1);
end;

################################################################
# find inverses in R
################################################################
ID := NewDictionary(false,true,A);
IC := NewInfoClass("TraceInverse");

ThinnedInverse := function(x)
  local a, b, c, d, t, u;

  if IsZero(Epsilon(x)) then return fail; fi;
  t := ThinnedReduce(x);
  if t in [one,one+A,one+B,one+C,one+D] then return t; fi;

  Info(IC,1,"Invert element with support ",Length(ThinnedSupport(t)));
### UGLY! we're probably stuck in some sort of recursion
  if Length(ThinnedSupport(t))>999 then Error("I'm probably stuck!"); fi;

  t := Thinned2Matrix(t);
  a := t[1][1]; b := t[1][2]; c := t[2][1]; d := t[2][2];

  if IsZero(Epsilon(a)) then
    t := ThinnedInverse(c);
    u := ThinnedInverse(b-a*t*d);
    t := [[-t*d*u,ThinnedInverse(c-d*ThinnedInverse(b)*a)],[u,-u*a*t]];
  else
    t := ThinnedInverse(d);
    u := ThinnedInverse(a-b*t*c);
    t := [[u,-u*b*t],[-t*c*u,ThinnedInverse(d-c*ThinnedInverse(a)*b)]];
  fi;
AddDictionary(ID,x,t);
  return Matrix2Thinned(t);
end;

CheckInverse := function(x)
  local z;
  z := ThinnedInverse(x);
  return [ThinnedIsZero(x*z-one),ThinnedIsZero(z*x-one)];
end;

################################################################
# multiplicative and nil order of element
################################################################
ThinnedOrder := function(x)
  local z;
  if IsZero(Epsilon(x)) then return fail;
  elif ThinnedIsZero(x-one) then return 1;
  else return 2*ThinnedOrder(ThinnedReduce(x^2));
  fi;
end;

ThinnedNilOrder := function(x)
  local z;
  if not IsZero(Epsilon(x)) then return infinity;
  elif ThinnedIsZero(x) then return 1;
  else return 2*ThinnedNilOrder(ThinnedReduce(x^2));
  fi;
end;

################################################################
# a basis of the nth power of the augmentation ideal,
# represented as a word in ABCD.
# W2Thinned returns the corresponding ring element.
################################################################
BasisOfDelta := function(n)
  local a, b, sigma, noa;

  sigma := w->Flat(List(w,function(s) if s='A' then return "ACA"; elif s='B' then return "D"; elif s='C' then return "B"; else return "C"; fi; end));
  noa := function(w)
    local i, j;
    if w[1]='A' then i := 2; else i := 1; fi;
    if w[Length(w)]='A' then j := Length(w)-1; else j := Length(w); fi;
    return w{[i..j]};
  end;
  if n=0 then return [""];
  elif n=1 then return ["A","B","D"];
  elif n=2 then return ["AB","BA","AD","DA"];
  elif n=3 then return ["ABA","ADA","BAB","BAD","DAB"];
  elif n=4 then return ["ABAB","ABAD","ADAB","BABA","BADA","DABA"];
  elif n=5 then return ["ABABA","ABADA","ADABA","BABAB","BABAD","BADAB","DABAB","DABAD"];
  elif IsOddInt(n) then
    a := List(BasisOfDelta((n-1)/2),w->noa(sigma(w)));
    b := List(BasisOfDelta((n+1)/2),w->noa(sigma(w)));
    return Union(List(a,x->Concatenation("A",x,"A")),b);
  else
    a := List(BasisOfDelta(n/2),w->noa(sigma(w)));
    return Union(List(a,x->Concatenation(x,"A")),List(a,x->Concatenation("A",x)));
  fi;
end;

W2Thinned := function(w)
  return Product(List(ShallowCopy(w),function(x)
    if x='A' then return A;
    elif x='B' then return B;
    elif x='C' then return C;
    elif x='D' then return D;
    fi; end));
end;

# inverse of 1+A+D+AB = (1+D)(1+AC)(1+ACAC)(1+A)
# inverse of 1+A+B+AD = (1+B)(1+AC)(1+ACAC)(1+A)

################################################################
# attempt to find inverse of element by long division
################################################################
ThinnedDegree := function(x)
  return Minimum(Union([infinity],List(ThinnedSupport(ThinnedReduce(x)),Length)));
end;

ThinnedLeadingTerm := function(x)
  local e, coeff, elm, deg;

  e := CoefficientsAndMagmaElements(ThinnedReduce(x));
  elm := e{[1,3..Length(e)-1]};
  coeff := e{[2,4..Length(e)]};
  deg := Minimum(Union([infinity],List(elm,Length)));
  e := Filtered([1..Length(elm)],i->Length(elm[i])=deg);

  return ElementOfMagmaRing(FamilyObj(x),ZeroCoefficient(x),coeff{e},elm{e});
end;

ThinnedLongDivisionInverse := function(x)
  local quo, rem, t;

  if ThinnedLeadingTerm(x)<>one then return fail; fi;

  quo := zero; rem := one; # running assumption: x*quo+rem = 1
  while rem<>zero do
    Info(IC,1,"Invert: remainder is ",rem);
    t := ThinnedLeadingTerm(rem);
    Info(IC,1,"Invert: leading term is ",t);
    quo := quo + t;
    rem := SuperReduce(rem + (x*t));
  od;
  return quo;
end;

Magma2Algebra := function(x)
  return ElementOfMagmaRing(FamilyObj(A),ZeroCoefficient(A),[One(F)],[x]);
end;

NthDegreeAlgebra := function(n)
  local b, sc, i, j, c, l;

  b := Union(List(Filtered(BasisOfDelta(n),w->w[1]='A'),x->ThinnedSupport(SuperReduce(W2Thinned(x)))));
  sc := EmptySCTable(Length(b),Zero(F));
  for i in Cartesian([1..Length(b)],[1..Length(b)]) do
    l := Thinned2Matrix(Magma2Algebra(b[i[1]]*b[i[2]]));
    c := CoefficientsAndMagmaElements(SuperReduce(l[1][1]+l[2][2]));
    l := [];
    for j in [1,3..Length(c)-1] do
      if Subword(c[j],1,1) in ThinnedSupport(A) then
        Add(l,One(F)); Add(l,Position(b,c[j]));
      fi;
    od;
    SetEntrySCTable(sc,i[1],i[2],l);
  od;
  return AlgebraByStructureConstants(F,sc);
end;

TwoOrder := function(x)
  local n, y;
  y := x;
  for n in [1..10] do
    if IsZero(y) then return n; fi;
    y := y*y;
  od;
  return infinity;
end;

################################################################
