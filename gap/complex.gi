#############################################################################
##
#W complex.gi                                               Laurent Bartholdi
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

# methods

if IsBound(COMPLEX_ROOTS_FR) then
    InstallMethod(RootsFloatOp, "for a list of coefficients",
            [IsList,IsPMComplex],
            function(l,tag)
        local r;
        r := COMPLEX_ROOTS_FR(List(l,z->NewFloat(IsPMComplex,z)));
        while r=fail do
            if Length(l)<=1 then return []; fi; # that's OK, no root
            Error("COMPLEX_ROOTS_FR returned Fail. Repent.");
        od;
        return List(r,z->Objectify(TYPE_PMCOMPLEX,z));
    end);
fi;

if IsBound(REAL_ROOTS_FR) then
    InstallMethod(RootsFloatOp, "for a list of coefficients",
            [IsList,IsIEEE754FloatRep],
            function(l,tag)
        local r;
        r := REAL_ROOTS_FR(List(l,z->NewFloat(IsIEEE754FloatRep,z)));
        while r=fail do
            if Length(l)<=1 then return []; fi; # that's OK, no root
            Error("REAL_ROOTS_FR returned Fail. Repent.");
        od;
        return r;
    end);
fi;

# printers

InstallMethod(String, [IsPMComplex], function(x)
    local s;
    s := ShallowCopy(String(x![1]));
    if x![2]>PMCOMPLEX_FAMILY!.rz then
        Add(s,'+');
    fi;
    if x![2]<>PMCOMPLEX_FAMILY!.rz then
        Append(s,String(x![2])); Add(s,'i');
    fi;
    Append(s,"_z");
    return s;
end);
InstallMethod(DisplayString, [IsPMComplex], function(x)
    local s;
    s := ShallowCopy(String(x![1]));
    if x![2]>PMCOMPLEX_FAMILY!.rz then
        Add(s,'+');
    fi;
    if x![2]<>PMCOMPLEX_FAMILY!.rz then
        Append(s,String(x![2])); Add(s,'i');
    fi;
    return s;
end);
InstallMethod(ViewString, [IsPMComplex], function(x)
    local i, s;
    s := ShallowCopy(ViewString(x![1]));
    i := x![2];
    if IsFinite(x![1]) then
        i := (i+x![1])-x![1]; # wipe out a very small imaginary part
    fi;
    if i>PMCOMPLEX_FAMILY!.rz then
        Add(s,'+');
    fi;
    if i<>PMCOMPLEX_FAMILY!.rz then
        Append(s,ViewString(i));
        Add(s,'i');
    fi;
    return s;
end);

# methods

InstallOtherMethod(IsZero, [IsPMComplex], x->IsZero(x![1]) and IsZero(x![2]));
InstallOtherMethod(IsOne, [IsPMComplex], x->IsOne(x![1]) and IsZero(x![2]));
InstallMethod(EQ, IsIdenticalObj, [IsPMComplex, IsPMComplex],
        function(x,y)
    return x![1]=y![1] and x![2]=y![2];
end);
InstallMethod(LT, "for complex numbers",
        [IsPMComplex,IsPMComplex],
        function(x,y)
    return x![1]<y![1] or (x![1]=y![1] and x![2]<y![2]);
end);

InstallMethod(RealPart, [IsPMComplex], x->x![1]);
InstallMethod(ImaginaryPart, [IsPMComplex], x->x![2]);
InstallMethod(ComplexConjugate, [IsPMComplex],
        x->Objectify(TYPE_PMCOMPLEX, [x![1],-x![2]]));
InstallMethod(Norm, [IsPMComplex], x->x![1]^2+x![2]^2);
InstallMethod(AbsoluteValue, [IsPMComplex], x->Sqrt(x![1]^2+x![2]^2));
InstallMethod(Argument, [IsPMComplex], x->ATAN2_MACFLOAT(x![2],x![1]));

InstallMethod(SUM, IsIdenticalObj, [IsPMComplex, IsPMComplex],
        function(x,y)
    return Objectify(TYPE_PMCOMPLEX, [x![1]+y![1], x![2]+y![2]]);
end);

InstallMethod(DIFF, IsIdenticalObj, [IsPMComplex, IsPMComplex],
        function(x,y)
    return Objectify(TYPE_PMCOMPLEX, [x![1]-y![1], x![2]-y![2]]);
end);

InstallMethod(AINV_MUT, [IsPMComplex], x->Objectify(TYPE_PMCOMPLEX, [-x![1],-x![2]]));

InstallOtherMethod(PROD, IsIdenticalObj, [IsPMComplex, IsPMComplex],
        function(x,y)
    return Objectify(TYPE_PMCOMPLEX, [x![1]*y![1]-x![2]*y![2],x![1]*y![2]+x![2]*y![1]]);
end);

InstallOtherMethod(INV, [IsPMComplex], function(x)
    local r;
    r := x![1]^2+x![2]^2;
    return Objectify(TYPE_PMCOMPLEX, [x![1]/r,-x![2]/r]);
end);

InstallOtherMethod(QUO, IsIdenticalObj, [IsPMComplex, IsPMComplex],
        function(x,y)
    return x*INV(y);
end);

InstallOtherMethod(POW, IsIdenticalObj, [IsPMComplex, IsPMComplex],
        function(x,y)
    local r, n, a;
    a := ATAN2_MACFLOAT(x![2],x![1]);
    n := Sqrt(x![1]^2+x![2]^2);
    r := n^y![1]*EXP_MACFLOAT(-y![2]*a);
    a := y![1]*a+y![2]*LOG_MACFLOAT(n);
    return Objectify(TYPE_PMCOMPLEX, [r*COS_MACFLOAT(a),r*SIN_MACFLOAT(a)]);
end);

InstallOtherMethod(POW, [IsPMComplex, IsScalar],
        function(x,y)
    local r, a;
    r := (x![1]^2+x![2]^2)^Float(y/2);
    a := ATAN2_MACFLOAT(x![2],x![1])*Float(y);
    return Objectify(TYPE_PMCOMPLEX, [r*COS_MACFLOAT(a),r*SIN_MACFLOAT(a)]);
end);

InstallOtherMethod(POW, [IsPMComplex, IsInt],
        function(x,n)
    local j, xpow, y;
    if n=0 then return NewFloat(IsPMComplex,1); elif n<0 then n := -n; x := INV(x); fi;
    if n>100 then TryNextMethod(); fi;
    y := NewFloat(IsPMComplex,1);
    while n<>0 do
        if IsOddInt(n) then y := y*x; fi;
        if n>1 then x := x*x; fi;
        n := QuoInt(n,2);
    od;
    return y;
end);

InstallMethod(Sqrt, [IsPMComplex],
        function(x)
    local r, a;
    r := Sqrt(Sqrt(x![1]^2+x![2]^2));
    a := ATAN2_MACFLOAT(x![2],x![1])*Float(1/2);
    return Objectify(TYPE_PMCOMPLEX, [r*COS_MACFLOAT(a),r*SIN_MACFLOAT(a)]);
end);

InstallMethod(Exp, [IsPMComplex],
        function(z)
    local r;
    r := EXP_MACFLOAT(z![1]);
    return Objectify(TYPE_PMCOMPLEX, [r*COS_MACFLOAT(z![2]),r*SIN_MACFLOAT(z![2])]);
end);

InstallOtherMethod(Random, [IsPMComplexCollection],
        function(D)
    if D=PMCOMPLEX_PSEUDOFIELD then
        return NewFloat(IsPMComplex,
                       Random(GlobalMersenneTwister, 0, 10^18)/10^18,
                       Random(GlobalMersenneTwister, 0, 10^18)/10^18);
    else
        TryNextMethod();
    fi;
end);

# complex.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
