#############################################################################
##
#W bisets.gi                                                Laurent Bartholdi
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  This file implements general bisets.
##
#############################################################################

BindGlobal("FRBISET_NEWTYPE@", function()
    return NewType(FRBISET_FAMILY,IsFRBiset and IsFRBisetByFRMachine);
end);

InstallMethod(BisetByFRMachine, [IsGroupFRMachine],
        function(M)
    local b, g;
    b := Objectify(FRBISET_NEWTYPE@(), rec(machine := M));
    g := SCGroup(M);
    SetLeftActingDomain(b,g);
    SetRightActingDomain(b,g);
    return b;
end);

InstallMethod(BisetByFRGroup, [IsFRGroup],
        function(G)
    local b;
    b := Objectify(FRBISET_NEWTYPE@(), rec(machine := UnderlyingFRMachine(G)));
    SetLeftActingDomain(b,G);
    SetRightActingDomain(b,G);
    return b;
end);


#E bisets.gi. . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
