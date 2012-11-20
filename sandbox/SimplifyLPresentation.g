
LoadPackage("kbmag");
LoadPackage("fr");

# This function tries to simplify an LPresentation
# - lpresgrp is an LPresentedGroup
# - maxeqns is a positive integer, which we give as an option the Knuth-Bendix
#   algorithm. The bigger maxeqns, the longer you will wait.
#   But your patience will (possibly) be rewarded by a yet simpler presentation
#   Try with maxeqns = 500 for example.
# - depth is a nonnegative integer which controls the number of endomorphisms to be considered.
#   Again, the bigger the depth, the longer you will wait.
#   Try with depth = 0 or depth = 1.
#   (we don't know any example where it was worth going further than depth = 1)

SimplifyLPresentation := function(lpresgrp, maxeqns, depth)
  local fgrp, frels, endos, irels, endosspheres, k, i, g, rws, reducedword;
  fgrp := FreeGroupOfLpGroup(lpresgrp);
  frels := ShallowCopy(FixedRelatorsOfLpGroup(lpresgrp));
  endos := EndomorphismsOfLpGroup(lpresgrp);
  irels := ShallowCopy(IteratedRelatorsOfLpGroup(lpresgrp));
  
  endosspheres := WordGrowth(Monoid(endos), rec(limit := depth, spheres := true));
  for k in [0..depth] do
    i := Length(irels);
    while i > 0 do
      g := fgrp / Concatenation(frels, Flat(ListX(irels{Concatenation([1..i-1],[i+1..Length(irels)])}, Flat(endosspheres{[1..k+1]}), \^)));
      rws := KBMAGRewritingSystem(g);
      OptionsRecordOfKBMAGRewritingSystem(rws).maxeqns := maxeqns;
      KnuthBendix(rws);
      reducedword := ReducedWord(rws, irels[i]);
      if reducedword = One(fgrp) then
        irels := irels{Concatenation([1..i-1],[i+1..Length(irels)])};
      else
        irels[i] := reducedword;
      fi;
      i := i - 1;
    od;
    i := Length(frels);
    while i > 0 do
      g := fgrp / Concatenation(frels{Concatenation([1..i-1],[i+1..Length(frels)])}, Flat(ListX(irels, Flat(endosspheres{[1..k+1]}), \^)));
      rws := KBMAGRewritingSystem(g);
      OptionsRecordOfKBMAGRewritingSystem(rws).maxeqns := maxeqns;
      KnuthBendix(rws);
      reducedword := ReducedWord(rws, frels[i]);
      if reducedword = One(fgrp) then
        frels := frels{Concatenation([1..i-1],[i+1..Length(frels)])};
      else
        frels[i] := reducedword;
      fi;
      i := i - 1;
    od;
  od;
  return LPresentedGroup(fgrp, frels, endos, irels);
end;
