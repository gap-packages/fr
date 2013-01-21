RequirePackage("fr");

SpectrumInfo := NewInfoClass("SpectrumInfo");

# SetInfoLevel(SpectrumInfo,1);

Spectrum := function(perms,nev,ncv)
    local p, i, j, n;
    
    n := LargestMovedPoint(perms);
    p := InputOutputLocalProcess(DirectoryCurrent(),"a.out",[]);
    WriteLine(p,Concatenation(" ",String(n)," ",String(nev)," ",String(ncv)," ",String(Size(perms)*n)," "));
    for i in [1..n] do for j in perms do
        WriteLine(p,Concatenation(" ",String(i)," ",String(i^j)));
    od; od;
    j := ReadLine(p);
    while Position(j,'\\')=fail do
        i := ReadLine(p);
        if i=fail then return [fail,j]; fi;
        if Position(i,'#')<>fail then Info(SpectrumInfo,1,NormalizedWhitespace(i));
        else Append(j,i); fi;
    od;
    j := [];
    for i in [1..nev] do
        Append(j,[FLOAT_STRING(ReadLine(p))]);
    od;
    CloseStream(p);
    return j;
end;

G := FRGroup("a=<b,b>(1,2)","b=<c,a>","c=<a,c>":MealyElement);
g := List([1..10],i->PermGroup(G,i));
