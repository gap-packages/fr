G := GrigorchukGroup;
Q := List([3..10],n->Image(EpimorphismPermGroupFrGroup(G,n)));

RW := function(G,n)
  local p, i;
  p := One(G);
  for i in [1..n] do p := p*Random(GeneratorsOfGroup(G)); od;
  return p;
end;

Test := function(G,len,num)
  local hits, i;
  hits := 0;
  for i in [1..num] do if RW(G,len)=One(G) then hits := hits+1; fi; od;
  return hits/num;
end;

Mat := List([1..4],d->List([5,10..100],len->Test(Q[d],len,5000)));

PrintTo("mat",Mat);
