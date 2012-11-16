DEBUG := rec(files := [],
             data := [],
             breakpoints := [],
             display := [],
             Read := function(name)
    local fileno, s, c, charno, lastcharno, lineno, lastlineno, semicolon;
    fileno := Position(DEBUG.files, name);
    if fileno=fail then
        Add(DEBUG.files, name);
        fileno := Length(DEBUG.files);
    fi;
    DEBUG.data[fileno] := ReadAll(InputTextFile(name));
    charno := 1;
    lastcharno := 1;
    lineno := 1;
    lastlineno := 1;
    s := "";
    semicolon := false;
    name := Filename(DirectoryTemporary(),"debug.g");
    for c in DEBUG.data[fileno] do
        Add(s,c);
        if c='\n' then
            lineno := lineno+1;
            if semicolon then
                AppendTo(name,Concatenation("DEBUG.BREAK(",String(fileno),",",String(lastlineno),",",String(lastcharno),",",String(charno),"); "));
                AppendTo(name,s);
                lastlineno := lineno;
                lastcharno := charno;
                s := "";
            fi;
        fi;
        charno := charno+1;
        semicolon := c=';';
    od;
    s := READ(s);
#    RemoveFile(name);
    return s;
end,
  BREAK := function(fileno,lineno,startpos,endpos)
    local i;
    if [fileno,lineno] in DEBUG.breakpoints then
        for i in DEBUG.display do
            Print(i,": ",EvalString(i),"\n");
        od;
        Error("Breakpoint in ",DEBUG.files[fileno],", line ",lineno,"\n",DEBUG.data[fileno]{[startpos..endpos]});
    fi;
end);

InstallMethod(Read, "with debugging support",
        [IsString],1000,
        function(name)
    local  readIndent, found;
    name := USER_HOME_EXPAND( name );
    if not IsReadableFile( name ) then
        Error( "file \"", name, "\" must exist and be readable" );
    fi;        
    readIndent := SHALLOW_COPY_OBJ( READ_INDENT );
    APPEND_LIST_INTR( READ_INDENT, "  " );
    InfoRead1( "#I", READ_INDENT, "Read( \"", name, "\" )\n" );
    if ValueOption("debug")=fail then
        found := READ(name);
    else
        found := DEBUG.Read(name);
    fi;
    READ_INDENT := readIndent;
    if found and READ_INDENT = ""  then
        InfoRead1( "#I  Read( \"", name, "\" ) done\n" );
    fi;
    if not found  then
        Error( "file \"", name, "\" must exist and be readable" );
    fi;
    return;
end);
