#############################################################################
##
#W  chapter-9-a.tst                FR Package               Laurent Bartholdi
##
#H  @(#)$Id$
##
#Y  Copyright (C) 2011,  Laurent Bartholdi
##
#############################################################################
##
##  This file tests the functions explained in chapter 9 of the manual,
##  that do not use the DLL.
##
#############################################################################

gap> START_TEST("fr:chapter 9 (1/2)");
gap> 
gap> Info(InfoFR,1,"9.2 Supporting rays");
#I  9.2 Supporting rays
gap> e := EquivalenceRelationPartition(ExternalAnglesRelation(2,5));
[ [ 1/31, 2/31 ], [ 1/15, 2/15 ], [ 3/31, 4/31 ], [ 1/7, 2/7 ],
  [ 5/31, 6/31 ], [ 1/5, 4/15 ], [ 7/31, 8/31 ], [ 9/31, 10/31 ],
  [ 1/3, 2/3 ], [ 11/31, 12/31 ], [ 2/5, 3/5 ], [ 13/31, 18/31 ],
  [ 3/7, 4/7 ], [ 14/31, 17/31 ], [ 7/15, 8/15 ], [ 15/31, 16/31 ],
  [ 19/31, 20/31 ], [ 21/31, 22/31 ], [ 5/7, 6/7 ], [ 11/15, 4/5 ],
  [ 23/31, 24/31 ], [ 25/31, 26/31 ], [ 13/15, 14/15 ], [ 27/31, 28/31 ],
  [ 29/31, 30/31 ] ]
gap> ForAll(e,p->SupportingRays(PolynomialIMGMachine(2,[p[1]]))[2][1][1]*2 in [p[1],p[2],p[1]+1,p[2]+1]);
true
gap> 
gap> Info(InfoFR,1,"Testing rabbit and z^2+i twists");
#I  Testing rabbit and z^2+i twists
gap> ri := PolynomialIMGMachine(2,[],[1/6]);
<FR machine with alphabet [ 1 .. 2 ] and adder FRElement(...,f4) on Group( [ f\
1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>
gap> model := StateSet(ri);
<free group on the generators [ f1, f2, f3, f4 ]>
gap> twist := GroupHomomorphismByImages(model,model,GeneratorsOfGroup(model),[model.1,model.2^(model.3*model.2),model.3^model.2,model.4]);
[ f1, f2, f3, f4 ] -> [ f1, f2^-1*f3^-1*f2*f3*f2, f2^-1*f3*f2, f4 ]
gap> 
gap> r := PolynomialIMGMachine(2,[1/7],[]);
<FR machine with alphabet [ 1 .. 2 ] and adder FRElement(...,f4) on Group( [ f\
1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>
gap> model := StateSet(r);
<free group on the generators [ f1, f2, f3, f4 ]>
gap> twist := GroupHomomorphismByImages(model,model,GeneratorsOfGroup(model),[model.1^(model.2*model.1),model.2^model.1,model.3,model.4]);
[ f1, f2, f3, f4 ] -> [ f1^-1*f2^-1*f1*f2*f1, f1^-1*f2*f1, f3, f4 ]
gap> rt := List([0..4],i->r*twist^i);
[ <FR machine with alphabet [ 1 .. 2 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <FR machine with alphabet [ 1 .. 2 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <FR machine with alphabet [ 1 .. 2 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <FR machine with alphabet [ 1 .. 2 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <FR machine with alphabet [ 1 .. 2 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]> ]
gap> m := PolynomialIMGMachine(3,[[3/4,1/12],[1/4,7/12]],[]);
<FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f3) on Group( [ f\
1, f2, f3 ] )/[ f3*f2*f1 ]>
gap> 
gap> Info(InfoFR,1,"Testing a folding");
#I  Testing a folding
gap> fold1 := NewIMGMachine("a=<,,b,,,B>(1,2,3)(4,5,6)","b=<,,b*a/b,,,B*A/B>","A=<,,b*a,,,B*A>(3,6)","B=(1,6,5,4,3,2)");
<FR machine with alphabet [ 1 .. 6 ] on Group( [ a, b, A, B ] )/[ a*B*A*b ]>
gap> 
gap> STOP_TEST( "chapter-9-a.tst", 10^10 );
fr:chapter 9 (1/2)
GAP4stones: 4681000

#E chapter-9-a.tst . . . . . . . . . . . . . . . . . . . . . . . . .ends here
