#############################################################################
##
#W read.g                                                   Laurent Bartholdi
##
#Y Copyright (C) 2006-2016, Laurent Bartholdi
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
ReadPackage("fr", "gap/perlist.gi");
ReadPackage("fr", "gap/frmachine.gi");
ReadPackage("fr", "gap/frelement.gi");
ReadPackage("fr", "gap/mealy.gi");
ReadPackage("fr", "gap/group.gi");
ReadPackage("fr", "gap/vhgroup.gi");
ReadPackage("fr", "gap/vector.gi");
ReadPackage("fr", "gap/linear.gi");
ReadPackage("fr", "gap/algebra.gi");
ReadPackage("fr", "gap/bisets.gi");
ReadPackage("fr", "gap/examples.gi");
ReadPackage("fr", "gap/cp.gi");
#############################################################################

# added to fix problems with loops
InstallMethod(\in,
         "(FR) default method, checking for <g> being among the generators",
          ReturnTrue,
         [IsFRElement, IsFRSemigroup], 1000,
   function ( g, G )
     if   g = One(G)
       or (IsFinite(GeneratorsOfGroup(G)) and g in GeneratorsOfGroup(G))
     then return true;
     else TryNextMethod(); fi;
end );

# added because ViewString cannot print ideals
Perform([[LeftActingRingOfIdeal,"left",GeneratorsOfLeftIdeal],
        [RightActingRingOfIdeal,"right",GeneratorsOfRightIdeal],
        [RightActingRingOfIdeal,"two-sided",GeneratorsOfTwoSidedIdeal]],
        function(data)
    InstallMethod(ViewString,
            [IsRing and Tester(data[1]) and Tester(data[3])],
            function(I)
        local s;
        s := "\>\><"; Append(s,data[2]); Append(s," ideal in \>\>");
        Append(s,ViewString(data[1](I)));
        Append(s,"\<,\< \>\>(");
        if HasDimension( I ) then
            Append(s,"dimension "); Append(s,String(Dimension(I)));
        else
            Append(s,String(Length(data[3](I)))); Append(s," generators");
        fi;
        Append(s,"\<\<\<\<)>");
        return s;
    end);
end);
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

DeclareAttribute("IsLevelTransitive", IsFRElement);
DeclareAttribute("IsLevelTransitive", IsFRGroup);
InstallMethod(IsLevelTransitive, [IsFRGroup], IsLevelTransitiveFRGroup);
InstallMethod(IsLevelTransitive, [IsFRElement], IsLevelTransitiveFRElement);
        
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
