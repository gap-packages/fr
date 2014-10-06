#############################################################################
##
##  This file tests the functions introduced by cp
##
#############################################################################
gap> START_TEST("frcp");
gap> n := InfoLevel(InfoFR);;
gap> SetInfoLevel(InfoFR,0);
gap> SetInfoLevel(InfoFRCP,1);
gap> Read(Filename(DirectoriesPackageLibrary("fr","tst"),"groups.g"));
gap> ForAll([1..Size(grps)],i->List(pairs[i],x->IsConjugate(grps[i],x[1],x[2])) = [true,false]);
true	
gap> List(sc_grps,x->List(pairs_sc,y->IsConjugate(x,y[1],y[2]))) = [[true,false,true],[true,false,true],[false,false,true],[false,false,false],[false,true,false],[false,true,false],[false,false,false],[false,false,false]];
true
## Test Method for branched groups
gap> L := GrigorchukConjugateBranchInit();
[ [ [ <Trivial Mealy element on alphabet [ 1 .. 2 ]>, a, 
          <Mealy element on alphabet [ 1 .. 2 ] with 6 states>, 
          <Mealy element on alphabet [ 1 .. 2 ] with 6 states> ], [  ], [  ], 
      [  ] ], 
  [ [  ], [ <Trivial Mealy element on alphabet [ 1 .. 2 ]>,,,, b,,,, c,,,, d ]
        , [  ], [  ] ], 
  [ [  ], [  ], 
      [ <Trivial Mealy element on alphabet [ 1 .. 2 ]>,,,, b,,,, c,,,, d ], 
      [  ] ], 
  [ [  ], [  ], [  ], 
      [ <Trivial Mealy element on alphabet [ 1 .. 2 ]>,,, 
          <Mealy element on alphabet [ 1 .. 2 ] with 6 states>, b,,, 
          <Mealy element on alphabet [ 1 .. 2 ] with 7 states>, c,,, 
          <Mealy element on alphabet [ 1 .. 2 ] with 7 states>, d,,, 
          <Mealy element on alphabet [ 1 .. 2 ] with 6 states> ] ] ]
gap> N:= GeneratorsOfGroup(GrigorchukGroup);
[ a, b, c, d ]
gap> InitConjugateForBranchGroups(GrigorchukGroup,N,L);
gap> IsConjugate(GrigorchukGroup,pairs[1][1][1],pairs[1][1][2]);
true
gap> SetInfoLevel(InfoFR,n);
gap> STOP_TEST( "cp.tst", 12*10^8 );
frcp
GAP4stones: 93000
