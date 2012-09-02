#############################################################################
##
#W p1.gd                                                    Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  This file declares code for P1 points
##
#############################################################################

## <#GAPDoc Label="P1Points">
## <ManSection>
##   <Filt Name="IsP1Point"/>
##   <Fam Name="P1PointsFamily"/>
##   <Func Name="P1Point" Arg="complex"/>
##   <Func Name="P1Point" Arg="real, imag" Label="ri"/>
##   <Func Name="P1Point" Arg="string" Label="s"/>
##   <Description>
##     P1 points are complex numbers or infinity;
##     fast methods are implemented to compute with them, and to apply
##     rational maps to them.
##     <P/>
##     The first filter recognizes these objects. Next, the family they
##     belong to. The next methods create a new P1 point.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="CleanedP1Point" Arg="p, prec"/>
##   <Returns><A>p</A>, rounded towards 0/1/infinity/reals at precision <A>prec</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Var Name="P1infinity"/>
##   <Var Name="P1one"/>
##   <Var Name="P1zero"/>
##   <Description>The south, north and 'east' poles of the Riemann sphere.</Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Antipode" Arg="p"/>
##   <Returns>The antipode of <A>p</A> on the Riemann sphere.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Barycentre" Arg="points ..."/>
##   <Returns>The barycentre of its arguments (which can also be a list of P1 points).</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Circumcentre" Arg="p, q, r"/>
##   <Returns>The centre of the smallest disk containing <A>p,q,r</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Distance" Arg="p, q"/>
##   <Returns>The spherical distance from <A>p</A> to <A>q</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Midpoint" Arg="p, q"/>
##   <Returns>The point between <A>p</A> to <A>q</A> (undefined if they are antipodes of each other).</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Sphere" Arg="v"/>
##   <Returns>The P1 point corresponding to <A>v</A> in <M>\mathbb R^3</M>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="SphereP1" Arg="p"/>
##   <Returns>The coordinates in <M>\mathbb R^3</M> of <A>p</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="SphereP1Y" Arg="p"/>
##   <Returns>The Y coordinate in <M>\mathbb R^3</M> of <A>p</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1XRatio" Arg="p q r s"/>
##   <Returns>The cross ratio of <A>p, q, r, s</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Filt Name="IsP1Map"/>
##   <Fam Name="P1MapsFamily"/>
##   <Description>
##     P1 maps are stored more efficiently than rational functions, but are
##     otherwise equivalent.
##     <P/>
##     The first filter recognizes these objects. Next, the family they
##     belong to. 
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="MoebiusMap" Arg="[sourcelist], destlist"/>
##   <Func Name="MoebiusMap" Arg="p, q, r, s, t, u" Label="6"/>
##   <Func Name="MoebiusMap" Arg="p, q, r" Label="3"/>
##   <Func Name="MoebiusMap" Arg="p, q" Label="2"/>
##   <Description>
##     These methods create a new P1 map. In the first case,
##     this is the Möbius transformation sending <A>p,q,r</A> to <A>P,Q,R</A>
##     respectively; in the second case, the map sending <A>p,q,r</A> to
##     <C>0,1,P1infinity</C> respectively; in the third case, the map sending
##     <A>p,q</A> to <C>0,P1infinity</C> respectively, of the form <M>(z-p)/(z-q)</M>.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Var Name="P1z"/>
##   <Description>The identity Möbius transformation.</Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="CleanedP1Map" Arg="map, prec"/>
##   <Returns><A>map</A>, with coefficients rounded using <A>prec</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="CoefficientsOfP1Map" Arg="map"/>
##   <Returns>Coefficients of numerator and denominator of <A>map</A>, lowest degree first.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapByCoefficients" Arg="numer, denom"/>
##   <Returns>The P1 map with numerator coefficients <A>numer</A> and denominator <A>denom</A>, lowest degree first.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1Path" Arg="p q"/>
##   <Returns>The P1 map sending <C>0</C> to <A>p</A> and <C>1</C> to <A>q</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="DegreeOfP1Map" Arg="map"/>
##   <Returns>The degree of <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1Image" Arg="map, p1point"/>
##   <Returns>The image of <A>p1point</A> under <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1PreImages" Arg="map, p1point"/>
##   <Returns>The preimages of <A>p1point</A> under <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapCriticalPoints" Arg="map"/>
##   <Returns>The critical points of <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapRational" Arg="rat"/>
##   <Returns>The P1 map given by the rational function <A>rat</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="RationalP1Map" Arg="map"/>
##   <Oper Name="RationalP1Map" Arg="indeterminate, map" Label="im"/>
##   <Returns>The rational function given by P1 map <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapSL2" Arg="mat"/>
##   <Returns>The Möbius P1 map given by the 2x2 matrix <A>mat</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="SL2P1Map" Arg="map"/>
##   <Returns>The matrix of the Möbius P1 map <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="SetP1Points" Arg="record [prec]"/>
##   <Description>
##     Installs a default implementation for P1 points. Fundamentally,
##     a P1 point is a complex number or infinity, with a few extra
##     methods. The argument <A>record</A> is the record describing
##     the floating-point implementation.
##     <P/>
##     Currently, one implementation (the default) is based on pairs
##     of IEEE754 floateans. It is fast, but is limited to 53 bits of
##     precision. It is loaded via <C>SetP1Points(PMCOMPLEX);</C>.
##     <P/>
##     Another implementation, in case the package
##     <Package>Float</Package> is available, is based on MPC complex
##     numbers. It offers unlimited precision, but is much slower. It is
##     loaded via <C>SetP1Points(MPC);</C> or <C>SetP1Points(MPC,prec);</C>.
##   </Description>
## </ManSection>
##
## <#/GAPDoc>

DeclareGlobalFunction("SetP1Points");

DeclareCategory("IsP1Point",IsObject);
DeclareCategory("IsIEEE754P1Point",IsP1Point);
BindGlobal("P1PointsFamily",NewFamily("P1PointsFamily",IsP1Point));
BindGlobal("TYPE_P1POINT",NewType(P1PointsFamily,IsP1Point and IsPositionalObjectRep));
BindGlobal("TYPE_IEEE754P1POINT",NewType(P1PointsFamily,IsIEEE754P1Point and IsDataObjectRep));

DeclareOperation("P1Point",[IsFloat]);
DeclareOperation("P1Point",[IsRat]);
DeclareOperation("P1Point",[IsInfinity]);
DeclareOperation("P1Point",[IsFloat,IsFloat]);
DeclareGlobalVariable("P1infinity");
DeclareOperation("P1INFINITY@",[IsP1Point]);
DeclareGlobalVariable("P1one");
DeclareGlobalVariable("P1zero");
DeclareOperation("P1Barycentre",[IsList]);
DeclareOperation("P1Barycentre",[IsP1Point]);
DeclareOperation("P1Barycentre",[IsP1Point,IsP1Point]);
DeclareOperation("P1Barycentre",[IsP1Point,IsP1Point,IsP1Point]);
DeclareAttribute("SphereP1",IsP1Point);
DeclareAttribute("SphereP1Y",IsP1Point);
DeclareAttribute("P1Sphere",IsList);
DeclareOperation("P1Distance",[IsP1Point,IsP1Point]);
DeclareOperation("P1Circumcentre",[IsP1Point,IsP1Point,IsP1Point]);
DeclareOperation("P1XRatio",[IsP1Point,IsP1Point,IsP1Point,IsP1Point]);
DeclareOperation("CleanedP1Point",[IsP1Point,IsFloat]);
DeclareOperation("P1Midpoint",[IsP1Point,IsP1Point]);
DeclareAttribute("P1Antipode",IsP1Point);
DeclareAttribute("P1Coordinate",IsP1Point);

################################################################
# p1 maps
################################################################
DeclareSynonym("IsP1Map",IsUnivariateRationalFunction and IsFloatRationalFunction);
DeclareCategory("IsIEEE754P1Map",IsP1Map);
BindGlobal("TYPE_IEEE754P1MAP", NewType(RationalFunctionsFamily(PMCOMPLEX_PSEUDOFIELD), IsIEEE754P1Map and IsDataObjectRep));

DeclareGlobalFunction("P1MapByCoefficients");
DeclareOperation("P1MAPBYCOEFFICIENTS2@",[IsObject,IsList,IsList]);
DeclareAttribute("CoefficientsOfP1Map",IsP1Map);
DeclareAttribute("AsP1Map",IsScalar);
DeclareOperation("P1MapSL2",[IsMatrix]);
DeclareAttribute("SL2P1Map",IsP1Map);
DeclareOperation("P1MapByZerosPoles",[IsList,IsList,IsP1Point,IsP1Point]);

DeclareOperation("MoebiusMap",[IsP1Point,IsP1Point]);
DeclareOperation("MoebiusMap",[IsP1Point,IsP1Point,IsP1Point]);
DeclareOperation("MoebiusMap",[IsP1Point,IsP1Point,IsP1Point,IsP1Point,IsP1Point,IsP1Point]);
DeclareOperation("MoebiusMap",[IsList]);
DeclareOperation("MoebiusMap",[IsList,IsList]);
DeclareOperation("P1Path",[IsP1Point,IsP1Point]);
DeclareGlobalVariable("P1z");
DeclareGlobalFunction("P1Monomial");
DeclareOperation("CleanedP1Map",[IsP1Map,IsFloat]);

DeclareAttribute("Primitive",IsP1Map);
#DeclareAttribute("Derivative",IsP1Map);
#DeclareAttribute("ComplexConjugate",IsP1Map);
DeclareAttribute("NumeratorP1Map",IsP1Map);
DeclareAttribute("DenominatorP1Map",IsP1Map);
DeclareSynonym("P1Image",ImageElm);
DeclareSynonym("P1PreImages",PreImagesElm);
DeclareAttribute("DegreeOfP1Map",IsP1Map);
DeclareSynonym("InverseP1Map",InverseGeneralMapping);
DeclareSynonym("CompositionP1Map",CompositionMapping2);
DeclareOperation("ImageElm",[IsP1Map,IsP1Point]);
DeclareOperation("PreImagesElm",[IsP1Map,IsP1Point]);
DeclareAttribute("CriticalPointsOfP1Map",IsP1Map);

DeclareOperation("P1INTERSECT@",[IsP1Map,IsP1Map,IsP1Map]);
DeclareOperation("P1ROTATION2@",[IsObject,IsList,IsObject]);

#############################################################################

#E p1.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
