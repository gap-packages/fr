#############################################################################
##
#W init.g                                                   Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2006, Laurent Bartholdi
##
#############################################################################
##
##  This file reads the declarations of the packages' new objects
##
#############################################################################

POSTHOOK@fr := []; # to be processed at the end

BindGlobal("@", rec()); # a record to store locals in the package

#############################################################################
##
#I Create info class to be able to debug loading
##
InfoFR := NewInfoClass("InfoFR");
SetInfoLevel(InfoFR, 1);
#############################################################################

#############################################################################
##
#R Read the declaration files.
##
ReadPackage("fr", "gap/helpers.gd");
ReadPackage("fr", "gap/complex.gd");
ReadPackage("fr", "gap/p1.gd");
ReadPackage("fr", "gap/perlist.gd");
ReadPackage("fr", "gap/trans.gd");
ReadPackage("fr", "gap/frmachine.gd");
ReadPackage("fr", "gap/frelement.gd");
ReadPackage("fr", "gap/mealy.gd");
ReadPackage("fr", "gap/group.gd");
ReadPackage("fr", "gap/vector.gd");
ReadPackage("fr", "gap/algebra.gd");
ReadPackage("fr", "gap/img.gd");
ReadPackage("fr", "gap/examples.gd");

ReadPackage("fr","hurwitz/gap/utils.gd");
ReadPackage("fr","hurwitz/gap/padicLift.gd");
ReadPackage("fr","hurwitz/gap/hurwitz.gd");

CallFuncList(function()
    local dirs, dll, w;
    dirs := DirectoriesPackagePrograms("fr");
    dll := Filename(dirs,"fr_dll.so");
    if dll=fail then
        dll := Filename(dirs[1],"fr_dll.so");
        for w in ["FIND_BARYCENTER","FIND_RATIONALFUNCTION",
                "C22P1POINT","P1POINT2C2","P1POINT2STRING","EQ_P1POINT",
                "P1SPHERE","SPHEREP1","SPHEREP1Y","P1BARYCENTRE",
                "P1ANTIPODE","P1MIDPOINT","P1DISTANCE","P1XRATIO",
                "CLEANEDP1POINT","P1CIRCUMCENTRE","LT_P1POINT",
                "P1MAPBYCOEFFICIENTS_IEEE754","P1INTERSECT_IEEE754",
                "MAT2P1MAP","P1MAP2MAT","P1MAP3","P1MAP2","P1PATH",
                "CLEANEDP1MAP","INVERSEP1MAP","COMPOSEP1MAP","P1IMAGE",
                "P1PREIMAGES","P1MAPCRITICALPOINTS","P1MAPCONJUGATE",
                "P1MAPBYZEROSPOLES","P1MAPPRIMITIVE","P1MAPDERIVATIVE",
                "P1MAPNUMER","P1ROTATION_IEEE754","P1MAPDENOM",
                "P1MAPISPOLYNOMIAL","STRINGS2P1POINT","DEGREEOFP1MAP",
                "P1MAPNUMERATOR","P1MAPDENOMINATOR","P1MAP_SUM","P1MAP_DIFF",
                "P1MAP_PROD","P1MAP_QUO","P1MAP_INV","P1MAP_AINV"] do
            CallFuncList(function(w)
                BindGlobal(w, function(arg)
                    Error("You need to compile ",dll," before using ",w,"\nYou may compile it with './configure && make' in ",PackageInfo("fr")[1].InstallationPath,"\n...");
                end);
            end,[w]);
        od;
        @.dll := false;
    else
        LoadDynamicModule(dll);
        @.dll := true;
    fi;
end,[]);

if not IsBound(IsLpGroup) then
    ForAll(["IsLpGroup","IsElementOfLpGroup","LPresentedGroup",
            "ElementOfLpGroup","SetEmbeddingOfAscendingSubgroup"], function(w)
        BIND_GLOBAL(w, fail);
        Add(POSTHOOK@fr,function() MAKE_READ_WRITE_GLOBAL(w); UNBIND_GLOBAL(w); end);
        return true;
    end);
fi;

InstallMethod(IsMatrixModule,[IsFRAlgebra],SUM_FLAGS,ReturnFalse);
# otherwise, bug causes SubmoduleNC(algebra,[]) to run indefinitely

#############################################################################

#E init.g . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
