# construct a word built from x,
# at given level, with structure m (in base 2)
#
Read("TreeGroups.g");
n := 7;
G := List([1..n],n->MakeGGroup(n));
NS := function(arg) return NormalClosure(G[n],Subgroup(G[n],arg)); end;
A := function(x) return Combine(2^(n-1),[(),Annihilate(x,2)]); end;
B := function(x) return Combine(2^(n-1),[Annihilate(x,2)^-1,Annihilate(x,2)]); end;
Wx := function(N,level,m)
  if m<0 or m>=2^level then Error("m must be between 0 and 2^level-1"); fi;
  if level=0 then
    return Comm(G[N].1,G[N].2);
  elif IsOddInt(m) then
    return Comm(G[N].1,Wx(N-1,level-1,QuoInt(m,2)));
  else
    return Wx(N-1,level-1,QuoInt(m,2));
  fi;
end;
W := function(m)
  if m=infinity then return ();
  else return Wx(n,LogInt(m,2),m-2^LogInt(m,2));
  fi;
end;
WW := function(m)
  return W(m)^2;
end;

#
# analyse the normal subgroup H.
#
Analyse := function(H)
  local j, k, l, K;
  j := First([1..2^(n-1)],n->WW(n) in H);
  for k in [1..2^(n-1)] do
    l := First(Combinations([1..j-1]),p->W(k)*Product(p,WW,()) in H);
    if l<>fail then
      K := NormalClosureInParent(Subgroup(G[n],[WW(j),W(k)*Product(l,WW,())]));
      if H=K then return [LogInt(IndexInParent(K),2),
	"W(k)*PROD W(l)^2,W(j)^2",k,l,j," rank ",LogInt(Index(H,CommutatorSubgroup(H,G[n])),2)]; fi;
      return [[LogInt(IndexInParent(H),2),LogInt(Index(H,K),2)],
	"W(k)*PROD W(l)^2,W(j)^2",k,l,j];
    fi;
  od;
end;
FastAnalyse := function(H)
  local j, k, l, K;
  j := First([1..2^(n-1)],n->WW(n) in H);
  for k in [1..2^(n-1)] do
    l := First(Combinations([1..j-1]),p->W(k)*Product(p,WW,()) in H);
    if l<>fail then return [LogInt(IndexInParent(H),2),"W(k)*PROD W(l)^2,W(j)^2",k,l,j]; fi;
  od;
end;

Binize := function(n)
  if n=1 then return "";
  elif IsOddInt(n) then return Concatenation("1",Binize((n-1)/2));
  else return Concatenation("0",Binize(n/2));
  fi;
end;
NameGroup := function(H)
  local l;
  l := FastAnalyse(H);
  if l[4]=[] then
    SetName(H,Concatenation(Binize(l[3]),"|",Binize(l[5])));
  else
    SetName(H,Concatenation(Binize(l[3]),Concatenation(List(l[4],i->Concatenation("*",Binize(i)))),"|",Binize(l[5])));
  fi;
end;

#
# construct subgroups
#
INITg := function()
  g := G[n];
  NS := function(arg) return NormalClosure(g,Subgroup(g,arg)); end;
  k := NS(W(1));
  k1 := NormalClosure(G[n-1],Subgroup(G[n-1],[Comm(G[n-1].1,G[n-1].2)]));
#  M := Set(Flat(List(Cartesian([1..10],[1..10]),x->NS(W(x[1]),W(x[2])^2))));
  I := function(l) return SortedList(List(l,x->Index(k,x))); end;
#  I5 := function() return SortedList(List(M5,x->Index(k5,x))); end;
end;

#
# lift subgroups from finite quotients
#
listsg := function(n)
  local pi;
  pi := NaturalHomomorphismByNormalSubgroup(k,n);
  return Filtered(List(NormalSubgroups(Image(pi)),x->PreImage(pi,x)),x->IsNormal(g,x));
end;

makeJK := function(j,k)
  return Filtered(List(Combinations([1..j-1]),function(l) local y;
	y := Analyse(NS(WW(j),W(k)*Product(l,WW,())));
	return [y[3]=j and y[4]=k and y[5]=l,y]; end),x->x[1]);
end;

#
# growth series
#
#for <a,b,c>:
#[1 3 5 8 13 21 31 46 69 102 151 222]
#for <a,b,d>:
#[1 3 5 8 12 17 25 37 54  79 116 170]
#for <a,c,d>:
#[1 3 5 8 12 17 25 37 53  75 107 152 213 295 407]
#for <a,b,c,d>, at level 5:
#[1 4 6 12 17 28 40 68 95 156
# 216 356 488 772 1054 1660 2218 3332 4450 6700
# 8773 12716 16538 23924 30161 40820 51444 69836 85796 112372
# 138708 182684 202815 216828 240300 267940 278748 282052 290753 293308
# 271007 228380 210627 186868 148625 95804 73495 53044 31765 13004
# 7305 4020 1495 368 193 96 19]

#
# fun with cohomology
#

F := GF(2);
Cohomology := function(g)
  local p, m;
  p := Range(IsomorphismPcGroup(g));
  m := GModuleByMats(List(Pcgs(p),x->IdentityMat(1,F)),F);
  return TwoCohomology(p,m).cohom;
end;

# for SylowSubgroup(SymmetricGroup(2^n),2), we get:
#n=1, 2, 3,  4,  5,  6
#d=1, 3, 7, 14, 25, 41 |guess: n^2
#
# for Grigorchuk group quotients, we get:
#n=1, 2, 3, 4,  5,  6
#d=1, 3, 7, 9, 11, 13
#
# for Supergroup quotients, we get:
#n=1, 2, 3,  4,  5,  6
#d=1, 3, 7, 14, 18, 22

# Cohomology(SylowSubgroup(SymmetricGroup(2^7),2));

#
# growth
#
# N := 10; series((1+X)^3*Product((1+X^(2*i))^2*(1+X^(2*i+1)),i=1..N),X,2*N);
