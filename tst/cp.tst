#############################################################################
##
##  This file tests the functions introduced by cp
##
#############################################################################
gap> START_TEST("frcp");
gap> n := InfoLevel(InfoFR);;
gap> sc_grps := [];;
gap> for alph in [[1,2],[1..3]] do
> 	for fil in [IsFRElement,IsFiniteStateFRElement,IsBoundedFRElement,IsFinitaryFRElement] do
> 	Add(sc_grps,FullSCGroup(alph,fil));
>    od;
> od;
gap> grps := [GrigorchukGroup,AleshinGroup,GuptaSidkiGroup];;

################################################################ init
# create pairs of group elements
gap> pairs := [];;
gap> for grp in grps do
> 	pairs_group:=[];
> 	n := Size(GeneratorsOfGroup(grp));
> 	Add(pairs_group,[grp.(2 mod n +1)^grp.(1 mod n +1),grp.(2 mod n +1)^grp.(3 mod n +1)]);
> 	Add(pairs_group,[grp.(1 mod n +1)^grp.(1 mod n +1),grp.(2 mod n +1)^grp.(3 mod n +1)]);
> 	Add(pairs,pairs_group);
> od;
gap> pairs_sc := [[AsGroupFRElement(AddingElement(2)),AsGroupFRElement(AddingElement(2))^-1],[AsGroupFRElement(AddingElement(3)),AsGroupFRElement(AddingElement(3))^-1],[AsGroupFRElement(AddingElement(2)),(AsGroupFRElement(AddingElement(2)))^FRElement([[[1],[AddingElement(2)]]],[()],[1])]];;

################################################################ tests
gap> ForAll([1..Size(grps)],i->List(pairs[i],x->IsConjugate(grps[i],x[1],x[2])) = [true,false]);
true
gap> List(sc_grps,x->List(pairs_sc,y->IsConjugate(x,y[1],y[2]))) = [[true,false,true],[true,false,true],[false,false,true],[false,false,false],[false,true,false],[false,true,false],[false,false,false],[false,false,false]];
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
#I  Degree: converting to Mealy element
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
>     RepSystem:=List(~.Branchstructure.group,x->PreImagesRepresentativeNC(~.Branchstructure.quo,x)))
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
