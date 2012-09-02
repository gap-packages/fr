gap> START_TEST("fr:p1");
gap> n := P1Point(0.0);
<.0e0>
gap> s := P1Antipode(n);
<∞>
gap> e := P1Point(1.0);;
gap> w := P1Point(-1.0);;
gap> f := P1Point(1.0i);;
gap> b := P1Point(-1.0i);;
gap> pts := [n,s,e,w,f,b,P1Point(1.e-3),P1Point(1.e-3i),P1Point(1.e-3+1.e-3i),P1Point(1.e-6i)];
[ <.0e0>, <∞>, <.1e1>, <-.1e1>, <.0e0+.1e1ⅈ>, <-.0e0-.1e1ⅈ>, <.1e-2>,
  <.0e0+.1e-2ⅈ>, <.1e-2+.1e-2ⅈ>, <.0e0+.1e-5ⅈ> ]
gap> P1Midpoint(s,n);
fail
gap> P1Midpoint(s,e);
<.241421e1>
gap> P1Midpoint(n,w);
<-.414214e0>
gap> P1Midpoint(pts[9],pts[10]);
<.5e-3+.5005e-3ⅈ>
gap> P1Midpoint(pts[10],pts[10]);
<.0e0+.1e-5ⅈ>
gap> List(pts,SphereP1);
[ [ .0e0, .0e0, .1e1 ], [ .0e0, .0e0, -.1e1 ], [ .1e1, .0e0, .0e0 ],
  [ -.1e1, -.0e0, .0e0 ], [ .0e0, .1e1, .0e0 ], [ -.0e0, -.1e1, .0e0 ],
  [ .2e-2, .0e0, .999998e0 ], [ .0e0, .2e-2, .999998e0 ],
  [ .2e-2, .2e-2, .999996e0 ], [ .0e0, .2e-5, .1e1 ] ]
gap> List(pts,x->P1Distance(x,P1Sphere(SphereP1(x))));
[ .0e0, .0e0, .0e0, .0e0, .0e0, .0e0, .0e0, .0e0, .0e0, .300927e-35 ]
gap> DelaunayTriangulation(pts);
<triangulation with 10 vertices, 48 edges and 16 faces>
gap> DelaunayTriangulation(pts,2.0);
<triangulation with 86 vertices, 504 edges and 168 faces>
gap> 

gap> z := Indeterminate(MPC_PSEUDOFIELD,"z");
z
gap> STOP_TEST("p1.tst", 10^8);
fr:p1
