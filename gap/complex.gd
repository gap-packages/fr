#############################################################################
##
#W complex.gd                                               Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  Poor man's complex numbers, implemented in GAP above IEEE754 machine-type
##
#############################################################################

## <#GAPDoc Label="complexnumbers">
## <ManSection>
##   <Filt Name="IsPMComplex"/>
##   <Fam Name="PMCOMPLEX_FAMILY"/>
##   <Var Name="PMCOMPLEX_PSEUDOFIELD"/>
##   <Var Name="PMCOMPLEX"/>
##   <Description>
##     A "poor man's" implementation of complex numbers, based on the
##     underlying 64-bit floating-point numbers in &GAP;.
##     <P/>
##     Strictly speaking, complex numbers do not form a field in &GAP;,
##     because associativity etc. do not hold. Still, a field is defined,
##     <C>PMCOMPLEX_FIELD</C>, making it possible to construct an indeterminate
##     and rational functions, to be passed to <Package>FR</Package>'s
##     routines.
##     <P/>
##     These complex numbers can be made the default floating-point numbers
##     via <C>SetFloats(PMCOMPLEX);</C>. They may then be entered as
##     standard floating-point numbers, with the suffix <C>_z</C>.
## <Example><![CDATA[
## gap> z := Indeterminate(PMCOMPLEX_FIELD,"z");
## z
## gap> (z+1/2)^5/(z-1/2);
## (z^5+2.5*z^4+2.5*z^3+1.25*z^2+0.3125*z+0.03125)/(z+(-0.5))
## gap> NewFloat(IsPMComplex,1,2);
## 1+2i
## gap> last^2;
## -3+4i
## gap> RealPart(last);
## -3
## gap> Norm(last2);
## 25
## gap> NewFloat(IsPMComplex,"1+2*I");
## 1+2i
## gap> RootsFloat(z^2-5);
## [ 2.23607, -2.23607 ]
## gap> RootsFloat(ListWithIdenticalEntries(80,1.0_z));
## [ 0.987688+0.156434i, 0.996917+0.0784591i, 0.996917-0.0784591i, 0.987688-0.156434i, 0.760406+0.649448i, 0.92388+0.382683i, 0.951057-0.309017i, 0.97237+0.233445i, 0.809017+0.587785i,
##   0.522499+0.85264i, 0.649448+0.760406i, 0.891007+0.45399i, 0.587785+0.809017i, 0.707107+0.707107i, 0.951057+0.309017i, 0.233445+0.97237i, 0.45399+0.891007i, 0.309017+0.951057i,
##   0.382683+0.92388i, 0.85264+0.522499i, -0.59719-0.608203i, -0.867574-0.11552i, -0.186972-0.990223i, -0.999006+0.318176i, -0.739308+0.0272973i, -0.432752-0.7287i, -0.672709+0.537561i,
##   0.156434+0.987688i, 0.295424-0.953359i, 0.588289-0.808509i, 0.455128-0.893999i, 0.0951213-1.01063i, 0.229628-0.939435i, -0.216054-0.95336i, -0.914152+0.49378i, 0.524052-0.853005i,
##   0.97237-0.233445i, -0.233486+0.972416i, 0.379514-0.92918i, 3.09131e-07+1.i, 0.182752-0.984684i, 0.891007-0.45399i, -0.0892207-1.01443i, 0.852641-0.522499i, 0.00247318-1.02032i,
##   0.92388-0.382683i, -0.585832+0.81608i, 0.809018-0.587792i, -0.656055+0.770506i, 0.760385-0.649467i, -0.452862+0.889692i, -0.0784562+0.996918i, 0.707015-0.707079i, 0.0784591+0.996917i,
##   -0.15643+0.987703i, -0.307608-0.969002i, 0.649377-0.760134i, -0.382904+0.92328i, -0.857704+0.573345i, -0.403754-0.946275i, -0.827986-0.648221i, -0.990655-0.396897i,
##   -0.929824-0.488558i, -0.671579-0.790133i, -0.886052-0.560249i, -1.05047-0.0873829i, -0.496236-0.900246i, -0.726008+0.713809i, -1.02514+0.223541i, -1.01722-0.277614i,
##   -0.585809-0.852796i, -0.518635+0.85364i, -1.04842+0.0255453i, -0.752485-0.724528i, -0.309225+0.951018i, -0.9612+0.409487i, -0.793651+0.646744i, -1.01735-0.194111i, -1.04161+0.124175i
##  ]
## gap> AsSortedList(List(last,AbsoluteValue));
## [ 0.739812, 0.847513, 0.852377, 0.861109, 0.875231, 0.967092, 0.977534, 0.998083, 0.998317, 0.998841, 0.99953, 0.999747, 0.999886, 0.999916, 0.999996, 1., 1., 1., 1., 1., 1., 1., 1.,
##   1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.00001, 1.00001, 1.00003, 1.00005, 1.00112, 1.0015, 1.00318, 1.0037, 1.00458, 1.00772, 1.01197,
##   1.01509, 1.01665, 1.01814, 1.01834, 1.02033, 1.0238, 1.02796, 1.02881, 1.03169, 1.03462, 1.0357, 1.03698, 1.03899, 1.04459, 1.04479, 1.04832, 1.04845, 1.04873, 1.04899, 1.04923,
##   1.05036, 1.05155, 1.0541, 1.05442, 1.0672 ]

## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>

DeclareCategory("IsPMComplex", IsFloat);
DeclareCategoryCollections("IsPMComplex");
DeclareCategoryCollections("IsPMComplexCollection");
BindGlobal("PMCOMPLEX_FAMILY", NewFamily("PMCOMPLEX_FAMILY", IsPMComplex));
SetIsUFDFamily(PMCOMPLEX_FAMILY,true);
PMCOMPLEX_FAMILY!.reals := IsIEEE754FloatRep;
PMCOMPLEX_FAMILY!.rz := NewFloat(PMCOMPLEX_FAMILY!.reals,0);
PMCOMPLEX_FAMILY!.ro := NewFloat(PMCOMPLEX_FAMILY!.reals,1);
BindGlobal("TYPE_PMCOMPLEX", NewType(PMCOMPLEX_FAMILY, IsPMComplex));

# global field, for polynomials

BindGlobal("PMCOMPLEX_PSEUDOFIELD",
        Objectify(NewType(CollectionsFamily(PMCOMPLEX_FAMILY),
                IsFloatPseudoField and IsAttributeStoringRep),rec()));
SetLeftActingDomain(PMCOMPLEX_PSEUDOFIELD,Rationals);
SetCharacteristic(PMCOMPLEX_PSEUDOFIELD,0);
# SetBaseField(PMCOMPLEX_PSEUDOFIELD,Rationals); # no such method seems to exist
SetDimension(PMCOMPLEX_PSEUDOFIELD,infinity);
SetSize(PMCOMPLEX_PSEUDOFIELD,infinity);
SetIsWholeFamily(PMCOMPLEX_PSEUDOFIELD,true);
SetName(PMCOMPLEX_PSEUDOFIELD,"PMCOMPLEX_PSEUDOFIELD");

# creators

INSTALLFLOATCREATOR("for pm complex", [IsPMComplex,IsScalar],
        function(filter,r)
    if IsPMComplex(r) then return r; fi;
    return Objectify(TYPE_PMCOMPLEX,[NewFloat(PMCOMPLEX_FAMILY!.reals,r),PMCOMPLEX_FAMILY!.rz]);
end);

InstallOtherMethod(NewFloat, "for pm complex", [IsPMComplex,IsScalar,IsScalar],
        function(filter,r,i)
    return Objectify(TYPE_PMCOMPLEX,[NewFloat(PMCOMPLEX_FAMILY!.reals,r),NewFloat(PMCOMPLEX_FAMILY!.reals,i)]);
end);

INSTALLFLOATCREATOR("for pm complex", [IsPMComplex,IsString],
        function(filter,s)
    local p, q, i, z;

    s := DifferenceLists(LowercaseString(s),Concatenation(WHITESPACE,"*"));
    if s in ["inf","infinity"] then return NewFloat(filter,infinity); fi;
    if s="nan" then return NewFloat(filter,PMCOMPLEX_FAMILY!.rz/PMCOMPLEX_FAMILY!.rz); fi;
    p := 1;     # start parsing a float here
    z := [PMCOMPLEX_FAMILY!.rz,PMCOMPLEX_FAMILY!.rz];
    while p <= Length(s) do
        i := 1;     # by default, real part
        q := p; # start parsing a float here
        if s[p] in "+-" then # sign
            p := p+1;
            if p>Length(s) then return fail; fi;
        fi;
        if s[p]='i' then
            i := 2; Remove(s,p); # imaginary part, zap it
            if p>Length(s) or s[p] in "+-" then Add(s,'1',p); fi;   # i+... = 1*i+...
        fi;
        while p<=Length(s) and s[p] in "0123456789." do p := p+1; od;
        if p<=Length(s) and s[p]='e' then # exponent
            p := p+1;
            if p<=Length(s) and s[p] in "+-" then p := p+1; fi;
            p := p+1;
            while p<=Length(s) and IsDigitChar(s[p]) do p := p+1; od;
        fi;
        if p<=Length(s) and s[p]='i' then
            if i=2 then return fail; fi; # two imaginaries
            i := 2; Remove(s,p);
        fi;
        if q>=p then return fail; fi; # no new characters
        q := NewFloat(PMCOMPLEX_FAMILY!.reals,s{[q..p-1]});
        if q=fail then return fail; fi; # something wrong
        z[i] := z[i] + q;
    od;
    return Objectify(TYPE_PMCOMPLEX,z);
end);

SetZero(PMCOMPLEX_FAMILY,NewFloat(IsPMComplex,0));
SetOne(PMCOMPLEX_FAMILY,NewFloat(IsPMComplex,1));
SetZero(PMCOMPLEX_PSEUDOFIELD,NewFloat(IsPMComplex,0));
SetOne(PMCOMPLEX_PSEUDOFIELD,NewFloat(IsPMComplex,1));

BindGlobal("PMCOMPLEX", rec(
    constants := rec(
        DIG := 15,
        VIEW_DIG := 6,
        MANT_DIG := 53,
        MAX_10_EXP := 308,
        MAX_EXP := 1024,
        MIN_10_EXP := -307,
        MIN_EXP := -1021,
        DECIMAL_DIG := 17,
        2IPI := NewFloat(IsPMComplex,0,2*ACOS_MACFLOAT(-1.0_l)),
        INFINITY := NewFloat(IsPMComplex,MACFLOAT_STRING("inf")),
        NAN := NewFloat(IsPMComplex,MACFLOAT_STRING("nan"))),
    filter := IsPMComplex,
    field := PMCOMPLEX_PSEUDOFIELD,
    reals := IEEE754FLOAT,
    creator := s->NewFloat(IsPMComplex,s),
    eager := 'z'));

SetFloats(PMCOMPLEX,false);

# complex.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
