CheckAH := function(f)
  local S,R,V,i,added;

  S := Source(f);
  R := Range(f);
  B := ShallowCopy(MappingGeneratorsImages(f)[1]);
  C := ShallowCopy(MappingGeneratorsImages(f)[2]);
  V := VectorSpace(LeftActingDomain(S),B);
  while V <> S do
    added := false;
    for i in Tuples([1..Length(B)],2) do
      if not B[i[1]]*B[i[2]] in V then
        Add(B,B[i[1]]*B[i[2]]);
        Add(C,C[i[1]]*C[i[2]]);
        V := VectorSpace(LeftActingDomain(S),B);
        added := true;
      fi;
    od;
    if IsRestrictedLieAlgebra(S) then for i in [1..Length(B)] do
      if not PthPowerImage(Basis(S),B[i]) in V then
        Add(B,PthPowerImage(Basis(S),B[i]));
        Add(C,PthPowerImage(Basis(R),C[i]));
        V := VectorSpace(LeftActingDomain(S),B);
        added := true;
      fi;
    od; fi;
    if not added then return fail; fi; # does not generate
  od;
  B := Basis(S,B);
  return [B,C];
  if IsRestrictedLieAlgebra(S) and not ForAll([1..Length(B)],
    i->LinearCombination(C,Coefficients(B,PthPowerImage(Basis(S),B[i])))=PthPowerImage(Basis(R),C[i])) then return false; fi;
  return ForAll(Tuples([1..Length(B)],2),
    i->LinearCombination(C,Coefficients(B,B[i[1]]*B[i[2]]))=C[i[1]]*C[i[2]]);
end;

G := MakeFrGroup("a=[b,c]","b=[c,b]","c=[a,a](1,2)");
q := List([1..15],i->PermGroupFrGroup(G,i));
r := DerivedSubgroup(q[10]);
k := function(n)
  return Group([[1,4],[0,1]]*ZmodnZObj(1,2^n),[[1,0],[4,1]]*ZmodnZObj(1,2^n),[[5,4],[-4,-3]]*ZmodnZObj(1,2^n));
end;
