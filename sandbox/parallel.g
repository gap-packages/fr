ls := function(M)
    local V, W, i, j;
    
    W := StateSet(M);
    repeat
        V := W;
        W := Subspace(V,[]);
        for i in Basis(V) do
            for j in M!.transitions do for j in j do
                W := ClosureLeftModule(W,i*j);
            od; od;
        od;
    until V=W;
    return VectorMachineNC(FamilyObj(M),MATRIX@FR(M!.transitions,x->List(Basis(W),b->Coefficients(Basis(W),b*x))),Basis(W)*M!.output);
end;

nuk := function(M)
    local l, oldl;
    
    M := Minimized(M);
    l := ls(M);
    repeat
        oldl := l;
        l := Minimized(ls(l*M));
    until l=oldl;
    return l;
end;

CONTROLLER := fail;
SLAVES := [];

BindGlobal("PAR_Stop", function()
    local r;
    
    IO_Pickle(CONTROLLER.write,"die");
    IO_kill(CONTROLLER.pid,9);
    IO_Close(CONTROLLER.read);
    IO_Close(CONTROLLER.write);
    CONTROLLER := fail;
end);

BindGlobal("PAR_Start", function(arg)
    local r, inputs, s;
    
    if CONTROLLER<>fail then
        Error("Controller already started, at pid ",CONTROLLER.pid);
    fi;
    r := [IO_pipe(),IO_pipe()];
    Add(r,IO_fork());
    if r[3]=0 then
        CONTROLLER := rec(pid := IO_getpid(),
                          readfd := r[2].toread,
                          read := IO_WrapFD(r[2].toread,false,false),
                          write := IO_WrapFD(r[1].towrite,false,false));
        # start slaves
        SLAVES := [];
        TODO := [];
        DONE := [];
        WAITID := fail;
        while true do
            inputs := Concatenation(List(SLAVES,r->r.readfd),[CONTROLLER.readfd]);
            r := IO_select(inputs,[],[],-1,-1);
            if r=0 or r=fail then continue; fi;
            if CONTROLLER.readfd in inputs then
                r := IO_Unpickle(CONTROLLER.read);
                if r="die" then
                    # kill slaves
                    IO_kill(CONTROLLER.pid,9);
                elif r="getdata" then
                    IO_Pickle(CONTROLLER.write,DONE);
                    DONE := [];
                elif r[1]="addslave" then
                    # add new slave
                elif r[1]="wait" then
                    WAITID := r{[2..3]};
                else
                    s := First(SLAVES,x->not x.idle);
                    Add(TODO,r);
                    # send to some slave, depending on priority
                fi;
            fi;
            for s in SLAVES do
                if s.readfd in inputs then
                    r := IO_Unpickle(s.read);
                    Add(DONE,r);
                    if r{[1..2]}=WAITID then # gap is waiting for a signal
                        IO_Pickle(CONTROLLER.write,true);
                        WAITID := fail;
                    fi;
                    if TODO=[] then
                        s.idle := true;
                    else
                        IO_Pickle(s.write,Remove(TODO)); # new task
                    fi;
                fi;
            od;
        od;
    else
        CONTROLLER := rec(pid := r[3],
                          readfd := r[1].toread,
                          read := IO_WrapFD(r[1].toread,false,false),
                          write := IO_WrapFD(r[2].towrite,false,false));
        InstallAtExit(PAR_Stop);
    fi;
end);

BindGlobal("PAR_Update", function()
    local d, obj;
    
    IO_Pickle(CONTROLLER.write,"getdata");
    for d in IO_Unpickle(CONTROLLER.read) do
        obj := DELAYED[d[1]];
        if IsList(obj[1]) then
            obj[1][d[2]] := d[3];
        else
            obj[1].(d[2]) := d[3];
        fi;
        obj[2] := obj[2]-1;
        if obj[2]=0 then
            Unbind(DELAYED[d[1]]);
        fi;
    od;
end);

BindGlobal("PAR_Push", function()
    # get taskid from list of delayed tasks; maybe add a line there
    IO_Pickle(CONTROLLER.write,[taskid,id,f,args,timeout]);
    return IO_Unpickle(CONTROLLER.write);
end);

BindGlobal("PAR_Wait", function()
    IO_Pickle(CONTROLLER.write,["wait",taskid,id]);
    IO_Unpickle(CONTROLLER.read); # blocking read
end);

DeclareCategory("IsDelayedList",IsList);
BindGlobal("DelayedListFamily",
        NewFamily("DelayedListsFamily",IsPeriodicList));
BindGlobal("TYPE_LIST_DELAYED",
        NewType(PeriodicListsFamily,IsPeriodicList));
DeclareOperation("DelayedList",[IsList,IsFunction]);

InstallMethod(DelayedList, [IsList,IsFunction],
        function(l,f)
    local i, m, pipe, pid;
    
    m := ValueOption("processes");
    if m=fail then m := 1; fi;
    pipe := [];
    pid := [];
    for i in [1..m] do
        Add(pipe,IO_pipe());
        Add(pid,IO_fork());
        if pid[i]=0 then
            # process to compute f(l(i))
            IO_kill(IO_getpid(),9);
        fi;
    od;
    return Objectify(TYPE_LIST_DELAYED,rec(data := [],
                   index := l, pipe := pipe, pid := pid));
end);

UPDATE_DELAYED_LIST := function(l)
    local i, r, s;
    while true do
        r := List(l!.pipe,x->x.toread);
        s := IO_select(r,[],[],0,0);
        if s=fail or s=0 then return; fi;
        for r in r do
            if r=fail then continue; fi;
            s := IO_Unpickle(IO_WrapFD(r,fail,fail));
            l!.data[s[1]] := s[2];
        od;
    od;
end;

InstallMethod(ViewObj, [IsDelayedList],
        function(l)
    UPDATE_DELAYED_LIST(l);
    Print("<delayed list, computed ",Number(l!.data),"/",Length(l!.index)," elements>\n");
end);

InstallMethod(\[\], [IsDelayedList,IsInt],
        function(l,i)
    UPDATE_DELAYED_LIST(l);
    if not IsBound(l!.data[i]) then
        # send request to l!.pipe
    fi;
    return l!.data[i];
end);
