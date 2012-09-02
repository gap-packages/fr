#############################################################################
##
#W  chapter-9-b.tst                 FR Package              Laurent Bartholdi
##
#H  @(#)$Id$
##
#Y  Copyright (C) 2011,  Laurent Bartholdi
##
#############################################################################
##
##  This file tests the functions explained in chapter 9 of the manual,
##  that use the DLL
##
#############################################################################

gap> START_TEST("fr:chapter 9 (2/2)");
gap> 
gap> Info(InfoFR,1,"12.8 P1 points");
#I  12.8 P1 points
gap> P1Distance(P1infinity,P1infinity);
0.
gap> 
gap> Info(InfoFR,1,"Shishikura-Tan Lei matings");
#I  Shishikura-Tan Lei matings
gap> SetFloats(IEEE754FLOAT);
gap> z := P1z;
<z>
gap> a := RootsFloat((z-1)*(3*z^2-2*z^3)+1);
[ 0.598631+0.565259i, -0.426536, 0.598631-0.565259i, 1.72927-2.22045e-16i ]
gap> c := RootsFloat((z^3+z)^3+z);
[ 0., 0.557573+0.540347i, -0.557573+0.540347i, -0.557573-0.540347i,
  0.557573-0.540347i, 0.264425+1.26049i, -0.264425+1.26049i,
  -0.264425-1.26049i, 0.264425-1.26049i ]
gap> am := List(a,a->IMGMachine((a-1)*(3*z^2-2*z^3)+1));
[ <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f3*f2*f1*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f2*f3*f1*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f3*f2*f1*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f3*f1*f2*f4 ]> ]
gap> cm := List(c,c->IMGMachine(z^3+c));
[ <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f2) on Group( [\
 f1, f2 ] )/[ f1*f2 ]>,
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f3*f1*f2*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f3*f1*f2*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f3*f1*f2*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f3*f1*f2*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f1*f3*f2*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f1*f3*f2*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f1*f3*f2*f4 ]>, 
  <FR machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Group( [\
 f1, f2, f3, f4 ] )/[ f1*f3*f2*f4 ]> ]
gap> m := ListX(am,cm,Mating);;
gap> RationalFunction(NewIMGMachine(am[2]));
<1.9353365831410638_z*z^3+(-6.5315742879608081_z)*z^2+5.2080077920476935_z*z+0\
.51788488659999687_z>
gap> RationalFunction(m[9+2]);
<((0.19853514848979054+0.094747862775783384i_z)*z^3+(0.36448217409192502+1.067\
1933890697718i_z)*z^2+(0.65828774719150207+0.56768707293267628i_z)*z+(-1.13054\
18372331943-1.8215142432570175i_z))/((0.42613189662774054+0.10853855675161406i\
_z)*z^3+(1.1298164009633644+1.950676369983551i_z)*z^2+(1.5115657119504351+0.85\
698452757134347i_z)*z+1._z)>
gap> 
gap> Info(InfoFR,1,"An obstructed mating");
#I  An obstructed mating
gap> RationalFunction(m[9+8]);
rec( machine := <FR machine with alphabet [ 1 .. 3 ] on Group(
    [ f1, f2, f3, g1, g2, g3 ] )/[ f2*f3*f1*g1*g3*g2 ]>,
  matrix := [ [ 1/2, 1 ], [ 1/2, 0 ] ], obstruction := [ f1*g1, f2^-1*g2^-1 ], 
  spider := <marked sphere on <triangulation with 9 vertices, 42 edges and 14 \
faces> marked by [ f1, f2, f3, g1, g2, g3 ] -> [ f1^-1*f4^-1, f3^-1*f2*f3, f3^\
-1*f5, f4*f1*f5^-1*f4^-1, f2^-1*f3, f4 ]> )
gap> 
gap> Info(InfoFR,1,"Testing Triangulations");
#I  Testing Triangulations
gap> if IsBound(MacFloat) then Float := MacFloat; fi;
gap> oct := List([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]],P1Sphere);;
gap> s := Sqrt(Float(1/3));;
gap> cube := List([[s,s,s],[s,s,-s],[s,-s,s],[-s,s,s],[s,-s,-s],[-s,s,-s],[-s,-s,s],[-s,-s,-s]],P1Sphere);;
gap> DelaunayTriangulation(cube);
<triangulation with 11 vertices, 54 edges and 18 faces>
gap> DelaunayTriangulation(cube{[1,5]});
<triangulation with 6 vertices, 24 edges and 8 faces>
gap> p := List([[0.,0.,1.],[0.,0.,-1.],SphereP1(P1Point(1.e-4)),SphereP1(P1Point(0.,1.e4))],P1Sphere);
[ <0+0i>, <P1infinity>, <0.0001+0i>, <-0+10000i> ]
gap> DelaunayTriangulation(p,100.);
<triangulation with 32 vertices, 180 edges and 60 faces>
gap> 
gap> Info(InfoFR,1,"Testing RationalFunction");
#I  Testing RationalFunction
gap> f := RationalFunction(PolynomialIMGMachine(2,[],[7/16]):param_unicritical);
<z^2+(-1.7712570233568821+0.066161509080687936i_z)>
gap> 
gap> Info(InfoFR,1,"Testing Pilgrim's obstructed blowup of the torus");
#I  Testing Pilgrim's obstructed blowup of the torus
gap> F := FreeGroup("a","b","c","d");
<free group on the generators [ a, b, c, d ]>
gap> Unbind(a); Unbind(b); Unbind(c); Unbind(d);
gap> AssignGeneratorVariables(F); o := One(F);;
#I  Assigned the global variables [ a, b, c, d ]
gap> M := FRMachine(F,[[c^-1,o,o,o,c],[o,o,o,d,d^-1],[a,o,o,a^-1,o],[b,o,d,a,c]],
>                   [(1,5)(2,4,3),(1,2)(4,5),(1,4)(2,3,5),()]);
<FR machine with alphabet [ 1 .. 5 ] on Group( [ a, b, c, d ] )>
gap> SetIMGRelator(M,d*c*b*a);
gap> RationalFunction(M);
rec(
  machine := <FR machine with alphabet [ 1 .. 5 ] on Group( [ a, b, c, d ] )/[\
 d*c*b*a ]>, matrix := [ [ 1 ] ], obstruction := [ a^-1*c^-1 ],
  spider := <marked sphere on <triangulation with 7 vertices, 30 edges and 10 \
faces> marked by [ a, b, c, d ] -> [ f1^-1*f2^-1, f3^-1*f1, f3, f2 ]> )
gap> 
gap> Info(InfoFR,1,"Testing mating of airplane with z^2+i");
#I  Testing mating of airplane with z^2+i
gap> m := Mating(PolynomialIMGMachine(2,[3/7],[]),PolynomialIMGMachine(2,[],[1/6]));
<FR machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2, f3, g1, g2, g3 ] )/[ \
f3*f2*f1*g3*g2*g1 ]>
gap> Unbind(f1); Unbind(f2); Unbind(f3); Unbind(g1); Unbind(g2); Unbind(g3); 
gap> AssignGeneratorVariables(StateSet(m));
#I  Assigned the global variables [ f1, f2, f3, g1, g2, g3 ]
gap> i := FreeGroup("f1","f2","f3","g1","x");
<free group on the generators [ f1, f2, f3, g1, x ]>
gap> tm := ChangeFRMachineBasis(m,[f1^-1*g2,One(StateSet(m))]);;
gap> inj := GroupHomomorphismByImages(i,StateSet(m),GeneratorsOfGroup(i),[f1^g2,f2,f3,g1,f1*g3/f1*g2]);
[ f1, f2, f3, g1, x ] -> [ g2^-1*f1*g2, f2, f3, g1, f1*g3*f1^-1*g2 ]
gap> m2 := SubFRMachine(tm,inj);
<FR machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2, f3, g1, x ] )/[ f3*f2\
*x*f1*g1 ]>
gap> RationalFunction(m2);
<((2.1173218773245877+0.092523780325534474i_z)*z^2+(3.6019898773482266+1.59353\
23696078603i_z)*z+(2.3547834193556909+1.5278264334816143i_z))/((-4.64791571561\
72919-5.8835287029321073i_z)*z+1._z)>
gap> 
gap> STOP_TEST( "chapter-9-b.tst", 10^10 );
fr:chapter 9 (2/2)
GAP4stones: 2498000

#E chapter-9-b.tst . . . . . . . . . . . . . . . . . . . . . . . . .ends here
