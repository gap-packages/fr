#############################################################################
##
#W p1.gi                                                    Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  This file contains helper code for P1 points, using generic floats
##
#############################################################################

InstallGlobalFunction(SetP1Points, function(arg)
    local floats, prec, eps;
    
    if Length(arg)=0 then
        floats := PMCOMPLEX;
        prec := floats.constants.MANT_DIG;
    else
        floats := arg[1];
        if Length(arg)=1 then
            if IsBound(floats.constants.MANT_DIG) then
                prec := floats.constants.MANT_DIG;
            else
                prec := 100; # default precision
            fi;
        else
            prec := arg[2];
        fi;
    fi;
    
    while not IsBound(floats.reals) do
        Error("The floats you provided don't seem to be complex numbers");
    od;
    
    SetFloats(floats.reals,prec,false);
    SetFloats(floats,prec,false);

    @.isc := floats.filter;
    @.isr := floats.reals.filter;
    @.field := floats.field;
    @.rfield := floats.reals.field;
    
    @.i := NewFloat(@.isc,0,1);
    @.pi := floats.reals.constants.PI;
    @.2ipi := floats.constants.2IPI;
    
    @.z := NewFloat(@.isc,0);
    @.o := NewFloat(@.isc,1);
    @.rz := NewFloat(@.isr,0);
    @.ro := NewFloat(@.isr,1);
    @.rinf := @.ro/@.rz;
    if IsBound(floats.reals.constants.EPSILON) then
        @.reps := floats.reals.constants.EPSILON;
    else
        eps := @.ro;
        repeat
            @.reps := eps;
            eps := eps / 2;
        until @.ro+eps = @.ro;
    fi;

    SetFloats(rec(eager := 'P', creator := s->P1Point(NewFloat(@.isc,s))),false);

    MakeReadWriteGlobal("P1zero");
    P1zero := P1Point(@.z);
    MakeReadOnlyGlobal("P1zero");
    
    MakeReadWriteGlobal("P1one");
    P1one := P1Point(@.o);
    MakeReadOnlyGlobal("P1one");
    
    MakeReadWriteGlobal("P1infinity");
    P1infinity := P1INFINITY@(P1zero);
    MakeReadOnlyGlobal("P1infinity");

    MakeReadWriteGlobal("P1z");
    P1z := P1MapByCoefficients([@.z,@.o],[@.o]);
    IsPolynomial(P1z);
    MakeReadOnlyGlobal("P1z");
    
    @.inum := IndeterminateNumberOfLaurentPolynomial(Indeterminate(@.field,"z"));
end);

################################################################
# p1 points, using package Floats
################################################################
InstallMethod(P1Point, "(FR) generic P1 point", [IsRat],
        z->P1Point(@.o*z));

InstallMethod(P1Point, "(FR) generic P1 point", [IsInfinity],
        z->P1Point(@.o,@.z));

InstallMethod(P1Point, "(FR) generic P1 point", [IsFloat],
        function(z)
    return Objectify(TYPE_P1POINT,[z]);
end);

InstallMethod(P1Point, "(FR) generic P1 point", [IsFloat,IsFloat],
        function(n,d)
    if IsZero(d) then
        return P1infinity;
    fi;
    return Objectify(TYPE_P1POINT,[n/d]);
end);

InstallMethod(P1INFINITY@, "(FR) generic P1 point", [IsP1Point],
        x->P1Point(@.o/@.z));

InstallMethod(P1Coordinate, "(FR) generic P1 point", [IsP1Point],
        x->x![1]);

InstallMethod(DisplayString, "(FR) generic P1 point", [IsP1Point],
        x->DisplayString(x![1]));

InstallMethod(ViewString, "(FR) generic P1 point", [IsP1Point],
        x->Concatenation("<",ViewString(x![1]),">"));

InstallMethod(String, "(FR) generic P1 point", [IsP1Point],
        x->Concatenation("P1Point(",String(x![1]),")"));

INSTALLPRINTERS@(IsP1Point);

InstallMethod(EQ, "(FR) generic P1 point", [IsP1Point,IsP1Point],
        function(p,q)
    return p![1]=q![1];
end);

InstallMethod(LT, "(FR) generic P1 point", [IsP1Point,IsP1Point],
        function(p,q)
    return p![1]<q![1];
end);

BindGlobal("C2_P1POINT@", function(p)
    p := p![1];
    if IsXInfinity(p) then
        return [@.o,@.z];
    elif Norm(p) <= @.ro then
        return [p,@.o];
    else
        return [@.o,Inverse(p)];
    fi;
end);

InstallMethod(P1Sphere, "(FR) generic coordinate list", [IsList],
        function(v)
    v := v / Sqrt(v*v);
    if v[3]>@.rz then
        return P1Point(NewFloat(@.isc,v[1],v[2]), v[3]+@.o);
    elif v[1]=@.rz and v[2]=@.rz then
        return P1infinity;
    else
        return P1Point(@.o-v[3], NewFloat(@.isc,v[1],-v[2]));
    fi;
end);

InstallMethod(SphereP1, "(FR) generic P1 point", [IsP1Point],
        function(p)
    local n;
    p := p![1];
    if IsXInfinity(p) then
        return [@.rz,@.rz,-@.ro];
    else
        n := Norm(p);
        return [2*RealPart(p),2*ImaginaryPart(p),@.ro-n] / (@.ro+n);
    fi;
end);

InstallMethod(SphereP1Y, "(FR) generic P1 point", [IsP1Point],
        function(p)
    p := p![1];
    if IsXInfinity(p) then
        return @.rz;
    else
        return 2*ImaginaryPart(p)/(@.ro+Norm(p));
    fi;
end);

InstallMethod(P1Antipode, "(FR) generic P1 point", [IsP1Point],
        function(p)
    p := p![1];
    if p=@.z then
        return P1infinity;
    elif IsXInfinity(p) then
        return P1zero;
    else
        return P1Point(-Inverse(ComplexConjugate(p)));
    fi;
end);

BindGlobal("CLEAN_CX@", function(z,prec)
    local r, i, changed;
    if IsXInfinity(z) then return z; fi;
    r := RealPart(z);
    i := ImaginaryPart(z);
    changed := false;
    if AbsoluteValue(i) < prec*AbsoluteValue(r) then
        if i<>@.rz then i := @.rz; changed := true; fi;
        if AbsoluteValue(r-@.ro) < prec and r<>@.ro then
            r := @.ro; changed := true;
        fi;
        if AbsoluteValue(r+@.ro) < prec and r<>-@.ro then
            r := -@.ro; changed := true;
        fi;
    fi;
    if AbsoluteValue(r) < prec*AbsoluteValue(i) then
        if r<>@.rz then r := @.rz; changed := true; fi;
        if AbsoluteValue(i-@.ro) < prec and i<>@.ro then
            i := @.ro; changed := true;
        fi;
        if AbsoluteValue(i+@.ro) < prec and i<>@.ro then
            i := -@.ro; changed := true;
        fi;
    fi;
    if changed then
        return NewFloat(@.isc,r,i);
    else
        return z;
    fi;
end);

InstallMethod(CleanedP1Point, "(FR) generic P1 point", [IsP1Point,IsFloat],
        function(p,prec)
    local z, w, n;
    z := p![1];
    w := CLEAN_CX@(z,prec);
    n := Norm(w);
    if n*2*prec*prec > @.ro then
        return P1infinity;
    elif n < 2*prec*prec then
        return P1one;
    elif IsIdenticalObj(z,w) then
        return p;
    else
        return P1Point(w);
    fi;
end);

InstallMethod(P1Barycentre, "(FR) generic list of P1 points", [IsList],
        list->P1Sphere(Sum(list,SphereP1)));

InstallMethod(P1Barycentre, "(FR) generic P1 point", [IsP1Point],
        p->p);

InstallMethod(P1Barycentre, "(FR) generic P1 point", [IsP1Point,IsP1Point],
        P1Midpoint);

InstallMethod(P1Barycentre, "(FR) generic P1 point", [IsP1Point,IsP1Point,IsP1Point],
        function(arg) return P1Barycentre(arg); end);
        
InstallMethod(P1Midpoint, "(FR) generic P1 point", [IsP1Point,IsP1Point],
        function(p,q)
    local a, d;
    p := p![1];
    q := q![1];
    if IsXInfinity(p) and IsXInfinity(q) then
        return P1infinity;
    elif IsXInfinity(p) then
        if q=@.z then return fail; fi;
        return P1Point(q*(@.o+Sqrt(@.o+@.o/Norm(q))));
    elif IsXInfinity(q) then
        if p=@.z then return fail; fi;
        return P1Point(p*(@.o+Sqrt(@.o+@.o/Norm(p))));
    fi;
    d := @.o + q*ComplexConjugate(p);
    if d=@.z then return fail; fi;
    a := Sqrt((@.o+Norm(p))/(@.o+Norm(q))*Norm(d));
    return P1Point(a*q+d*p,a+d);
end);

InstallMethod(P1Distance, "(FR) generic P1 point", [IsP1Point,IsP1Point],
        function(p,q)
    local v, d;
    p := p![1];
    q := q![1];
    if p=q then
        return @.rz;
    elif IsXInfinity(p) then
        v := @.ro/AbsoluteValue(q);
    elif IsXInfinity(q) then
        v := @.ro/AbsoluteValue(p);
    else
        d := @.o + q*ComplexConjugate(p);
        if d=@.z then
            v := @.rinf;
        else
            v := AbsoluteValue((p-q)/d);
        fi;
    fi;
    return 2*Atan(v);
end);

InstallMethod(P1XRatio, "(FR) generic P1 point", [IsP1Point,IsP1Point,IsP1Point,IsP1Point],
        function(p1,p2,p3,p4)
    p1 := C2_P1POINT@(p1);
    p2 := C2_P1POINT@(p2);
    p3 := C2_P1POINT@(p3);
    p4 := C2_P1POINT@(p4);
    return (p1[1]*p3[2]-p3[1]*p1[2])
           / (p2[1]*p3[2]-p3[1]*p2[2])
           * (p2[1]*p4[2]-p4[1]*p2[2])
           / (p1[1]*p4[2]-p4[1]*p1[2]);
end);

InstallMethod(P1Circumcentre, "(FR) generic P1 point", [IsP1Point,IsP1Point,IsP1Point],
        function(a,b,c)
    local p, q, v, i, d, centre;
    v := [C2_P1POINT@(a),C2_P1POINT@(b),C2_P1POINT@(c)];
    p := @.z;
    q := @.z;
    for i in [1..3] do
        a := v[i]; b := v[i mod 3 +1]; c := v[(i+1) mod 3 +1];
        p := p + Norm(a[1])*b[2]*c[2]*ComplexConjugate(b[1]*c[2]-c[1]*b[2]);
        q := q + a[1]*b[2]*ComplexConjugate(b[1]*a[2])*(Norm(c[1])+Norm(c[2]));
    od;
    q := (q - ComplexConjugate(q)) / 2;
    
    if p=@.z then
        centre := @.z;
    elif Norm(p)<Norm(q) then
        centre := -ComplexConjugate(p) / (q + Sqrt(q*q-Norm(p)));
    else
        centre := (-q + Sqrt(q*q-Norm(p))) / p;
    fi;
    
    d := AbsoluteValue(centre*v[1][2] - v[1][1]) / AbsoluteValue(ComplexConjugate(centre)*v[1][1] + v[1][2]);
    
    if d > @.ro then
        d := @.ro/d;
        centre := -@.o/ComplexConjugate(centre);
    fi;
    
    return [P1Point(centre), 2*Atan(d)];
end);

###############################################################################
### P1 maps
###############################################################################
InstallGlobalFunction(P1MapByCoefficients, function(arg)
    local map, numer, denom;
    while not Length(arg) in [1..2] or not ForAll(arg,IsHomogeneousList) do
        Error("Argument should be one or two lists of coefficients, not ",arg);
    od;
    numer := arg[1];
    if Length(arg)=2 then
        denom := arg[2];
    else
        denom := [One(numer[1])];
    fi;
    map := P1MAPBYCOEFFICIENTS2@(numer[1],numer,denom);
    return map;
end);

InstallMethod(AsP1Map, [IsP1Map], function(rat)
    if IsIdenticalObj(TypeObj(rat),TypeObj(P1z)) then
        return rat;
    else
        return CallFuncList(P1MapByCoefficients,List([NumeratorOfRationalFunction(rat),DenominatorOfRationalFunction(rat)],CoefficientsOfUnivariatePolynomial));
    fi;
end);

InstallMethod(P1MAPBYCOEFFICIENTS2@, "(FR) generic P1 map", [IsFloat,IsList,IsList],
        function(dummy,numer,denom)
    return UnivariateRationalFunctionByCoefficients(@.field,numer*@.o,denom*@.o,0,1);
end);

InstallMethod(CoefficientsOfP1Map, "(FR) generic P1 map", [IsP1Map],
        function(map)
    local c;
    map := CoefficientsOfUnivariateRationalFunction(map);
    c := [Concatenation(ListWithIdenticalEntries(map[3],@.z),map[1]),
          Concatenation(ListWithIdenticalEntries(-map[3],@.z),map[2])];
    c[3] := Maximum(Length(c[1]),Length(c[2]))-1;
    while Length(c[1])<=c[3] do Add(c[1],@.z); od;
    while Length(c[2])<=c[3] do Add(c[2],@.z); od;
    return c;
end);

InstallMethod(MoebiusMap, "(FR) for generic images of 0,1,infinity",
        [IsP1Point,IsP1Point,IsP1Point],
        function(p,q,r)
    local pq, qr;
    p := C2_P1POINT@(p);
    q := C2_P1POINT@(q);
    r := C2_P1POINT@(r);
    pq := q[1]*p[2]-p[1]*q[2];
    qr := r[1]*q[2]-q[1]*r[2];
    return P1MapByCoefficients([p[1]*qr,r[1]*pq],[p[2]*qr,r[2]*pq]);
end);

InstallMethod(MoebiusMap, "(FR) for generic images of 0,infinity",
        [IsP1Point,IsP1Point],
        function(p,q)
    p := C2_P1POINT@(p);
    q := C2_P1POINT@(q);
    return P1MapByCoefficients([p[1],q[1]],[p[2],q[2]]);
end);

InstallMethod(P1Path, "(FR) Möbius transformation 0->p, 1->q, infty->P1Antipode(p)",
        [IsP1Point,IsP1Point],
        function(p,q)
    local r, pq, qr;
    p := C2_P1POINT@(p);
    q := C2_P1POINT@(q);
    r := [-ComplexConjugate(p[2]),ComplexConjugate(p[1])];
    pq := q[1]*p[2]-p[1]*q[2];
    qr := r[1]*q[2]-q[1]*r[2];
    return P1MapByCoefficients([p[1]*qr,r[1]*pq],[p[2]*qr,r[2]*pq]);
end);

InstallMethod(MoebiusMap, "(FR) for 3 points in a list",
        [IsHomogeneousList],
        l->CallFuncList(MoebiusMap,l));

InstallMethod(MoebiusMap, "(FR) for 3 source points and 3 target points",
        [IsHomogeneousList,IsHomogeneousList],
        function(src,dst)
    return CompositionP1Map(CallFuncList(MoebiusMap,dst),
                   InverseP1Map(CallFuncList(MoebiusMap,src)));
end);

InstallMethod(MoebiusMap, "(FR) for 3 source points and 3 target points",
        [IsP1Point,IsP1Point,IsP1Point,IsP1Point,IsP1Point,IsP1Point],
        function(a,b,c,d,e,f)
    return CompositionP1Map(MoebiusMap(d,e,f),InverseP1Map(MoebiusMap(a,b,c)));
end);

InstallMethod(P1MapSL2, "(FR) for a matrix", [IsMatrix],
        function(mat)
    while not Length(mat)=2 and ForAll(mat,x->Length(x)=2) do
        Error("Argument ",mat," should be a 2x2 matrix");
    od;
    return P1MapByCoefficients(mat[1]{[2,1]},mat[2]{[2,1]});
end);

InstallMethod(SL2P1Map, "(FR) generic P1 map", [IsP1Map],
        function(map)
    while DegreeOfP1Map(map)<>1 do
        Error("Argument ",map," should be a degree-1 P1 map");
    od;
    return CoefficientsOfP1Map(map){[1,2]}{[2,1]};
end);

InstallGlobalFunction(P1Monomial, function(d)
    local m;
    m := NullMat(2,AbsoluteValue(d)+1);
    if d>=0 then
        m[1][d+1] := 1; m[2][1] := 1;
    else
        m[1][1] := 1; m[2][-d+1] := 1;
    fi;
    m := P1MapByCoefficients(m[1],m[2]);
    if d>=0 then IsPolynomial(m); else IsLaurentPolynomial(m); fi; # force true
    return m;
end);

InstallMethod(CleanedP1Map, "(FR) generic P1 map", [IsP1Map,IsFloat],
        function(map,prec)
    local norm, maxnorm, i, j, deg;
    
    map := List(CoefficientsOfP1Map(map),ShallowCopy);
    deg := Length(map[1]);
    for i in [1..2] do
        norm := List(map[i],Norm);
        maxnorm := prec*Maximum(norm);
        for j in [1..deg] do
            if norm[j] < maxnorm then map[i][j] := @.z; fi;
        od;
    od;
    i := 1; while IsZero(map[2][i]) do i := i+1; od;
    maxnorm := map[2][i];
    for i in [1..2] do
        for j in [1..deg] do
            map[i][j] := CLEAN_CX@(map[i][j]/maxnorm,prec);
        od;
    od;
    return P1MapByCoefficients(map[1],map[2]);
end);

InstallOtherMethod(InverseP1Map, "(FR) generic P1 map", [IsP1Map],
        function(map)
    local c;
    while DegreeOfP1Map(map)<>1 do
        Error("I don't know how to invert a map of degree >1");
    od;
    c := CoefficientsOfP1Map(map);
    return P1MapByCoefficients([-c[1][1],c[2][1]],[c[1][2],-c[2][2]]);
end);

InstallMethod(ConjugatedP1Map, "(FR) generic P1 map", [IsP1Map,IsP1Map],
        function(map,mobius)
    return CompositionP1Map(InverseP1Map(mobius),map,mobius);
end);

BindGlobal("POLY_DER@", function(coeff)
    local v, i;
    v := [];
    for i in [1..Length(coeff)-1] do
        Add(v,i*coeff[i+1]);
    od;
    return v;
end);

BindGlobal("POLY_MUL@", function(coeffa, coeffb)
    local v, i, j, sum, dega, degb;
    v := [];
    dega := Length(coeffa)-1;
    degb := Length(coeffb)-1;
    for i in [0..dega+degb] do
        sum := @.z;
        for j in [0..i] do
            if j <= dega and i-j <= degb then
                sum := sum + coeffa[j+1]*coeffb[i-j+1];
            fi;
        od;
        Add(v,sum);
    od;
    return v;
end);

BindGlobal("POLY_MULZMINUSB@", function(coeff, b) # multiply by (z-b)
    local v, i;
    if IsXInfinity(b) then
        v := ShallowCopy(coeff);
        Add(v,@.z);
    else
        v := [-b*coeff[1]];
        for i in [1..Length(coeff)-1] do
            Add(v,coeff[i] - b*coeff[i+1]);
        od;
        Add(v,coeff[Length(coeff)]);
    fi;
    return v;
end);

BindGlobal("POLY_EVAL@", function(coeff, x) # evaluate coeff at x
    local v, i, deg;
    deg := Length(coeff)-1;
    v := coeff[deg+1];
    for i in [deg,deg-1..1] do
        v := v*x + coeff[i];
    od;
    return v;
end);

BindGlobal("YLOP_EVAL@", function(coeff, x) # evaluate x^deg*coeff at 1/x
    local v, i;
    v := coeff[1];
    for i in [2..Length(coeff)] do
        v := v*x + coeff[i];
    od;
    return v;
end);

BindGlobal("DPOLY_EVAL@", function(coeff, x) # evaluate coeff' at x
    local v, i, deg;
    deg := Length(coeff)-1;
    v := deg*coeff[deg+1];
    for i in [deg-1,deg-2..1] do
        v := v*x + i*coeff[i+1];
    od;
    return v;
end);

InstallOtherMethod(CompositionP1Map, [IsP1Map,IsP1Map],
        function(map2,map1)
    local pow1, i, j, deg, num, den;
    map1 := CoefficientsOfP1Map(map1);
    map2 := CoefficientsOfP1Map(map2);
    deg := map1[3]*map2[3];
    
    # pow1[i] := numer1^i denom1^(deg2-i)
    pow1 := [map1[2],map1[1]];
    for i in [2..map2[3]] do
        for j in [i,i-1..1] do
            pow1[j+1] := POLY_MUL@(pow1[j],map1[1]);
        od;
        pow1[1] := POLY_MUL@(pow1[1],map1[2]);
    od;
    return P1MapByCoefficients(map2[1]*pow1,map2[2]*pow1);
end);

InstallOtherMethod(CompositionP1Map, [IsP1Map, IsP1Map, IsP1Map],
        function(map3,map2,map1)
    return CompositionP1Map(map3,CompositionP1Map(map2,map1));
end);

BindGlobal("P1MAP_EVAL@", function(num,den,x)
    if Norm(x) <= @.ro then
        return POLY_EVAL@(num,x) / POLY_EVAL@(den,x);
    else
        x := @.o / x;
        return YLOP_EVAL@(num,x) / YLOP_EVAL@(den,x);
    fi;
end);

InstallMethod(P1Image, "(FR) generic P1 map", [IsP1Map,IsP1Point],
        function(map,z)
    map := CoefficientsOfP1Map(map);
    return P1Point(P1MAP_EVAL@(map[1],map[2],z![1]));
end);

InstallOtherMethod(POW, "(FR) generic P1 map", [IsP1Point,IsP1Map],
        function(z,map)
    return P1Image(map,z);
end);

InstallOtherMethod(CallFuncList, "(FR) generic P1 map", [IsP1Map,IsList],
        function(map,zs)
    zs := List(zs,z->P1Image(map,z));
    if Length(zs)=1 then
        return zs[1];
    else
        return zs;
    fi;
end);

InstallMethod(P1PreImages, "(FR) generic P1 map", [IsP1Map,IsP1Point],
        function(map,z)
    local roots;
    map := CoefficientsOfP1Map(map);
    z := C2_P1POINT@(z);
    roots := List(RootsFloat(map[1]*z[2]-map[2]*z[1]),P1Point);
    while Length(roots)<map[3] do
        Add(roots,P1infinity);
    od;
    return roots;
end);

InstallMethod(CriticalPointsOfP1Map, "(FR) generic P1 map", [IsP1Map],
        function(map)
    local roots;
    map := CoefficientsOfP1Map(map);
    roots := List(RootsFloat(POLY_MUL@(POLY_DER@(map[1]),map[2])-POLY_MUL@(POLY_DER@(map[2]),map[1])),P1Point);
    while Length(roots)<2*map[3]-2 do
        Add(roots,P1infinity);
    od;
    return roots;
end);

InstallMethod(P1MapByZerosPoles, "(FR) generic P1 map", [IsList,IsList,IsP1Point,IsP1Point],
        function(zeros,poles,src,dst)
    # construct a rational map with specified zeros and poles, and sending
    # src to dst
    local den, num, z;
    
    while Length(zeros)<>Length(poles) or not ForAll(zeros,IsP1Point) or not ForAll(poles,IsP1Point) do
        Error("P1MapByZerosPoles: first 2 arguments must be lists of same length of zeros and poles");
    od;
    
    num := [@.o]; for z in zeros do num := POLY_MULZMINUSB@(num,z![1]); od;
    den := [@.o]; for z in poles do den := POLY_MULZMINUSB@(den,z![1]); od;
    
    if dst=P1infinity then
        num := num * src![1];
    else
        num := num * dst![1] / P1MAP_EVAL@(num,den,src![1])![1];
    fi;

    return P1MapByCoefficients(num,den);
end);

InstallMethod(P1INTERSECT@, "(FR) generic P1 map", [IsP1Map,IsP1Map,IsP1Map],
        function(gamma,ratmap,delta)
    # compute the (t,u) in [t0,1]x[0,1] such that gamma(t) = ratmap(delta(u)).
    # returns a list of [t,u,Im(gamma^-1*ratmap*delta)'(u),gamma(t),delta(u)]
    # gamma, delta are Möbius transformations, and ratmap is a rational map.
    local poly, roots, intersect, eps, t, u, z, tu;
#MARKTIME@(1);
    ratmap := CoefficientsOfP1Map(CompositionP1Map(InverseP1Map(gamma),ratmap,delta)); # 1.5ms
#MARKTIME@(2);    
    poly := POLY_MUL@(ratmap[1],List(ratmap[2],ComplexConjugate)); # 350mus
#MARKTIME@(3);
    roots := RootsFloat(List(poly,ImaginaryPart)); # 60ms -- 45ms for complex, 45mus for native, 80mus for complex
#MARKTIME@(4);    
    eps := NewFloat(@.isr,10^-8);
    
    intersect := [];
    for u in roots do
        if not @.isr(u) or u < -eps or u > 1+eps then continue; fi;
        z := P1MAP_EVAL@(ratmap[1],ratmap[2],u);
        # t = gamma^-1*ratmap*delta(u)
        t := RealPart(z);
        if IsXInfinity(z) or ImaginaryPart(z) < -@.ro or ImaginaryPart(z) > @.ro
           or t < -eps or t > 1+eps then
            # in fact, ImaginaryPart(z) is microscopic; just avoid infinity
            continue;
        fi;
        tu := [t,u,,P1Image(gamma,P1Point(z)),P1Image(delta,P1Point(u))];
        z := DPOLY_EVAL@(poly,u);
        z := ImaginaryPart(z) / AbsoluteValue(z); # direction of approach
        if z < -eps then
            tu[3] := -1;
        elif z > eps then
            tu[3] := 1;
        else
            tu[3] := 0;
        fi;
        Add(intersect,tu);
    od;
#MARKTIME@(5);
    return intersect;
end);

InstallMethod(P1ROTATION2@, "(FR) generic P1 points", [IsP1Point,IsList,IsObject],
        function(dummy,points,extra)
    # find a Möbius transformation that sends the last of points to
    # P1infinity, and either
    # - matches points and extra as well as possible, if extra is a list;
    # - does a dilatation around infinity of amplitude extra, if extra is real,
    #   and is only a rotation, otherwise.
    local i, p, moeb, proj, oldproj, theta, n;
    
    p := points[Length(points)]![1];
    if IsXInfinity(p) then
        moeb := [[@.z,@.o],[@.o,@.z]];
    elif Norm(p) <= @.ro then
        moeb := [[@.o,ComplexConjugate(p)],[-p,@.o]];
    else
        p := @.o/p;
        moeb := [[ComplexConjugate(p),@.o],[@.o,-p]];
    fi;

    proj := [];
    for p in points do
        p := P1MAP_EVAL@(moeb[1],moeb[2], p![1]);
        if IsXInfinity(p) then
            Add(proj,@.z);
        else
            Add(proj,2*p / (@.o + Norm(p)));
        fi;
    od;
    
    theta := @.z;
    
    if IsList(extra) then
        n := @.z;
        for i in [1..Length(points)] do
            p := extra[i]![1];
            if IsXInfinity(p) then
                oldproj := @.z;
            else
                oldproj := 2*p / (@.o + Norm(p));
            fi;
            theta := theta + ComplexConjugate(proj[i])*oldproj;
            n := n + Norm(proj[i]);
        od;
        theta := theta / n;
        if n=@.z or Norm(theta) < 7/10*@.ro then # no good rotation
            theta := @.z;
        fi;
    fi;
    
    if theta=@.z then
        # hard... as last resort, just force the point of largest
        # projection to be on the positive real axis
        p := @.ro/10;
        theta := @.o;
        for i in [1..Length(points)] do
            n := Norm(proj[i]);
            if n > p then p := n; theta := ComplexConjugate(proj[i]); fi;
        od;
    fi;
    theta := theta / AbsoluteValue(theta); # make it of norm 1
    
    if IsFloat(extra) then
        theta := theta * extra;
    fi;
    
    return P1MapByCoefficients(theta*moeb[1],moeb[2]);
end);

InstallMethod(NumeratorP1Map, "(FR) generic P1 map", [IsP1Map],
        NumeratorOfRationalFunction);

InstallMethod(DenominatorP1Map, "(FR) generic P1 map", [IsP1Map],
        DenominatorOfRationalFunction);

InstallMethod(DegreeOfP1Map, "(FR) for a rational function", [IsP1Map],
        f->Maximum(DegreeOfUnivariateLaurentPolynomial(
                NumeratorOfRationalFunction(f)),
                DegreeOfUnivariateLaurentPolynomial(
                        DenominatorOfRationalFunction(f))));

InstallOtherMethod(ComplexConjugate, "(FR) for a univariate rational function",
        [IsP1Map],
        function(f)
    local c;
    c := CoefficientsOfUnivariateRationalFunction(f);
    return UnivariateRationalFunctionByExtRepNC(FamilyObj(f),List(c[1],ComplexConjugate),List(c[2],ComplexConjugate),c[3],IndeterminateNumberOfUnivariateRationalFunction(f));
end);

InstallMethod(Primitive, "(FR) for a univariate polynomial",
        [IsP1Map and IsLaurentPolynomial],
        function(f)
    local d, i, c;

    d := CoefficientsOfLaurentPolynomial(f);
    if d[1]=[] then # easy case: primitive of 0-Polynomial
        return f;
    fi;
    c := [];
    for i in [1..Length(d[1])]  do
        if i=-d[2] then
            if not IsZero(d[1][i]) then TryNextMethod(); fi; # has log(x) term
            c[i] := d[1][i];
        else
            c[i] := d[1][i]/(i+d[2]);
        fi;
    od;
    return LaurentPolynomialByCoefficients(CoefficientsFamily(FamilyObj(f)),c,
                   d[2]+1,IndeterminateNumberOfUnivariateRationalFunction(f));
end);

###########################################################################

#E p1_mpc.gi . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
