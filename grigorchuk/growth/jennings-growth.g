RequirePackage("treegp");

G := Image(EpimorphismPermGroupFrGroup(GrigorchukGroup,7));
#L := JenningsSeries(G);
L := LowerCentralSeries(G);
n := Length(L)-1;
deg := function(x)
  if IsOne(x) then return n+1; else return InTower(x,L); fi;
end;
epi := List([1..n],i->NaturalHomomorphismByNormalSubgroup(L[i],L[i+1]));
quo := List(epi,f->Filtered(Image(f),x->not IsOne(x)));
rep := List([1..n],i->function(x) return Representative(PreImages(epi[i],x)); end);
quo[n+1] := [()];
wt := List(quo,x->List(x,x->999999));
wt[1] := List(quo[1],x->1);
set_wt := function(x,w)
  local d;
  d := deg(x);
  wt[d][Position(quo[d],x^epi[d])] := w;
end;
get_wt := function(x)
  local d;
  if IsOne(x) then return 0; fi;
  d := deg(x);
  return wt[d][Position(quo[d],x^epi[d])];
end;
iter := function()
  local i, j, k, l, x;
  for i in [1..n] do for j in [i..n] do
    for k in [1..Length(quo[i])] do for l in [1..Length(quo[j])] do
      x := Comm(rep[i](quo[i][k]),rep[j](quo[j][l]));
      if get_wt(x) > wt[i][k]+wt[j][l] then set_wt(x,wt[i][k]+wt[j][l]); fi;
    od; od;
  od; od;
end;

a := G.1; b := G.2; c := G.3; d := G.4;
A := x->Comm(x,a); B := x->Comm(x,b); C := x->Comm(x,c);
D := x->Comm(x,d); BC := x->Comm(x,b)*Comm(x,c);
BC := x->Comm(x,b)*Comm(x,c); BD := x->Comm(x,b)*Comm(x,d);
CD := x->Comm(x,c)*Comm(x,d); BCD := x->Comm(x,b)*Comm(x,c)*Comm(x,d);
ops := [A,B,C,D,BC,BD,CD,BCD];
x := function(l)
  local g, i;
  g := a;
  for i in l do g := ops[i](g); od;
  return [g,deg(g)];
end;

ab := B(a); ac := C(a); ad := D(a);             #2
aw := W(a); aba := A(ab);                       #3
awa := A(aw);                                   #4
awb := B(aw); awc := C(aw); awd := D(aw);       #5
awba := A(awb); awca := A(awc); awda := A(awd); #6



aww := W(aw);
