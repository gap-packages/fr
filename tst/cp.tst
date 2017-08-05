#############################################################################
##
##  This file tests the functions introduced by cp
##
#############################################################################
gap> START_TEST("frcp");
gap> n := InfoLevel(InfoFR);;
gap> SetInfoLevel(InfoFR,0);
gap> SetInfoLevel(InfoFRCP,1);
gap> Read(Filename(DirectoriesPackageLibrary("fr","tst"),"groups.g"));
gap> ForAll([1..Size(grps)],i->List(pairs[i],x->IsConjugate(grps[i],x[1],x[2])) = [true,false]);
true
gap> List(sc_grps,x->List(pairs_sc,y->IsConjugate(x,y[1],y[2]))) = [[true,false,true],[true,false,true],[false,false,true],[false,false,false],[false,true,false],[false,true,false],[false,false,false],[false,false,false]];
true

## Test Method for branched groups
gap> SetIsBranched(GrigorchukGroup,true);
gap> IsConjugate(GrigorchukGroup,pairs[1][1][1],pairs[1][1][2]);
#I  Init FRBranchGroupConjugacyData
#I  Finished Init FRBranchGroupConjugacyData
true

# Setting the recursion end manually. 
gap> SetFRBranchGroupConjugacyData(GuptaSidkiGroup,
>  rec(initial_conj_dic:=NewDictionary([One(GuptaSidkiGroup),One(GuptaSidkiGroup)],true),
>     Branchstructure:=BranchStructure(GuptaSidkiGroup),
>     RepSystem:=List(~.Branchstructure.group,x->PreImagesRepresentative(~.Branchstructure.quo,x)))
>  );
gap> CallFuncList(function(a,t) 
>             local G,D,g,h;
>             G:= GuptaSidkiGroup;
>             D:= FRBranchGroupConjugacyData(G).initial_conj_dic;
>             for g in [a,a^2,t,t^2] do
>               for h in [a,a^2,t,t^2] do
>                if g<>h then
>                 AddDictionary(D,[g,h],[]);
>                fi;
>               od;
>             od;
>             AddDictionary(D,[a,a],[One(G),a,a^2]);
>             AddDictionary(D,[a^2,a^2],[One(G),a,a^2]);
>             AddDictionary(D,[t,t],[One(G),,,t,,,t^2]);
>             AddDictionary(D,[t^2,t^2],[One(G),,,t,,,t^2]);
>           end,GeneratorsOfGroup(GuptaSidkiGroup)
>   );
gap> SetIsBranched(GuptaSidkiGroup,true);
gap> IsConjugate(GuptaSidkiGroup,GuptaSidkiGroup.1^GuptaSidkiGroup.2,GuptaSidkiGroup.1);
true

#
gap> SetInfoLevel(InfoFR,n);
gap> STOP_TEST( "cp.tst", 15*10^8 );
