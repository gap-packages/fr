#############################################################################
##
#W read.g                                                   Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2006, Laurent Bartholdi
##
#############################################################################
##
##  This file reads the implementations, and in principle could be reloaded
##  during a GAP session.
#############################################################################

#############################################################################
##
#R Read the install files.
##
ReadPackage("fr", "gap/helpers.gi");
ReadPackage("fr", "gap/complex.gi");
ReadPackage("fr", "gap/p1.gi");
ReadPackage("fr", "gap/p1_ieee754.gi");
if @.dll then
    SetP1Points(PMCOMPLEX);
elif IsPackageMarkedForLoading("float","") then
    if IsBound(MPC) then
        SetP1Points(MPC,100);
    elif IsBound(CXSC) then
        SetP1Points(CXSC);
    else
        Info(InfoPackageLoading, 1, "You installed the Float package but compiled neither the Float nor the FR DLL. That's probably not what you wanted. I'll disable the P1Points code.");
        @.ro := 1.0_l; # minimal defaults to shut up warnings -- won't be usable.
        @.o := NewFloat(IsPMComplex,1);
        @.reps := IEEE754FLOAT.constants.EPSILON;
    fi;
else
    Info(InfoPackageLoading, 2, "You didn't install the Float package, and didn't compile the FR DLL. I'll disable the P1Points code.");
    @.ro := 1.0_l; # minimal defaults to shut up warnings -- won't be usable.
    @.o := NewFloat(IsPMComplex,1);
    @.reps := IEEE754FLOAT.constants.EPSILON;
fi;
ReadPackage("fr", "gap/perlist.gi");
ReadPackage("fr", "gap/trans.gi");
ReadPackage("fr", "gap/frmachine.gi");
ReadPackage("fr", "gap/frelement.gi");
ReadPackage("fr", "gap/mealy.gi");
ReadPackage("fr", "gap/group.gi");
ReadPackage("fr", "gap/vhgroup.gi");
ReadPackage("fr", "gap/vector.gi");
ReadPackage("fr", "gap/linear.gi");
ReadPackage("fr", "gap/algebra.gi");
ReadPackage("fr", "gap/img.gi");
ReadPackage("fr", "gap/examples.gi");
#############################################################################

#############################################################################
##
#X install shortcuts
##
INSTALL@ := function()
    CallFuncList(function(arg)
        local s;
        for s in arg do
            if IsBoundGlobal(s) then
                Info(InfoFR,2,Concatenation("Removing cybersquatter `",s,"'"));
                if IsReadOnlyGlobal(s) then MakeReadWriteGlobal(s); fi;
                UnbindGlobal(s);
            fi;
        od;
    end, ["Nucleus","Decomposition"]);

    DeclareOperation("Nucleus", [IsFRSemigroup]);
    InstallMethod(Nucleus, [IsFRSemigroup], NucleusOfFRSemigroup);
    DeclareOperation("Nucleus", [IsFRMachine]);
    InstallMethod(Nucleus, "(FR) for an FR machine", [IsFRMachine], NucleusOfFRMachine);

    DeclareOperation("Decomposition", [IsFRMachine]);
    InstallMethod(Decomposition, "(FR) for an FR element", [IsFRElement], DecompositionOfFRElement);
end;
#############################################################################

DeclareAttribute("Alphabet", IsFRObject);
DeclareAttribute("Alphabet", IsFRSemigroup);
DeclareAttribute("Alphabet", IsFRAlgebra);
InstallMethod(Alphabet, [IsFRObject], AlphabetOfFRObject);
InstallMethod(Alphabet, [IsFRSemigroup], AlphabetOfFRSemigroup);
InstallMethod(Alphabet, [IsFRAlgebra], AlphabetOfFRAlgebra);

if IsBound(Nucleus) and FLAG2_FILTER(Nucleus)=0 then
    DeclareOperation("Nucleus", [IsFRMachine]);
    DeclareOperation("Nucleus", [IsFRSemigroup]);
else
    DeclareAttribute("Nucleus", IsFRMachine);
    DeclareAttribute("Nucleus", IsFRSemigroup);
fi;
InstallMethod(Nucleus, [IsFRMachine], NucleusOfFRMachine);
InstallMethod(Nucleus, [IsFRSemigroup], NucleusOfFRSemigroup);

if IsBound(Decomposition) and FLAG2_FILTER(Decomposition)=0 then
    DeclareOperation("Decomposition", [IsFRElement]);
else
    DeclareAttribute("Decomposition", IsFRElement);
fi;
InstallMethod(Decomposition, [IsFRElement], DecompositionOfFRElement);

while not IsEmpty(POSTHOOK@fr) do Remove(POSTHOOK@fr)(); od;
Unbind(POSTHOOK@fr);

if IsBound(IO_Pickle) then
    ReadPackage("fr","gap/pickle.g");
else
    if not IsBound(IO_PkgThingsToRead) then
        IO_PkgThingsToRead := [];
    fi;
    Add(IO_PkgThingsToRead, ["fr","gap/pickle.g"]);
fi;

#E read.g . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
