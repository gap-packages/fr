#############################################################################
##
#W linear-gbnp.gi                                           Laurent Bartholdi
##
##
#Y Copyright (C) 2007, Laurent Bartholdi
##
#############################################################################
##
##  Methods to be installed only if GBNP is present
##
#############################################################################

InstallMethod(FRMachineRWS, "(FR) for an algebra machine",
        [IsLinearFRMachine and IsAlgebraFRMachineRep],
        function(M)
    local rws;
    rws := rec(free := M!.free, gbasis := [], gbasiscopy := []);

    rws.restart := function()
        rws.gbasis := ShallowCopy(rws.gbasiscopy);
    end;

    rws.commit := function()
        rws.gbasiscopy := ShallowCopy(rws.gbasis);
    end;

    rws.reduce := function(x)
        return NP2GP(StrongNormalFormNP(GP2NP(x),rws.gbasis),M!.free);
    end;
    #!!! maybe work purely in the GBNP format, to speed up?

    rws.addrule := function(x)
        Add(rws.gbasis,GP2NP(x));
        rws.gbasis := GBNP.ReducePol(rws.gbasis);
    end;

    return rws;
end);

BindGlobal("ALGEBRAISZERO@", function(M,x)
    local rws, todo, i;

    rws := NewFRMachineRWS(M);
    todo := NewFIFO([x]);
    for x in todo do
        x := rws.reduce(x);
        if not IsZero(x) then
            if not IsZero(SUBS@(x,M!.output)) then return false; fi;
            rws.addrule(x);
            x := SUBS@(x,M!.transitions);
            for i in x do Append(todo,i); od;
        fi;
    od;
    rws.commit();
    return true;
end);

