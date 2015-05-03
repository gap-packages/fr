sc_grps := [];

for alph in [[1,2],[1..3]] do
	for fil in [IsFRElement,IsFiniteStateFRElement,IsBoundedFRElement,IsFinitaryFRElement] do
		Add(sc_grps,FullSCGroup(alph,fil));
	od;
od;

grps := [GrigorchukGroup,AleshinGroup,GuptaSidkiGroup];
#################
#pairs of group elements
#################
pairs := [];
for grp in grps do
	pairs_group:=[];
	n := Size(GeneratorsOfGroup(grp));
	Add(pairs_group,[grp.(2 mod n +1)^grp.(1 mod n +1),grp.(2 mod n +1)^grp.(3 mod n +1)]);
	Add(pairs_group,[grp.(1 mod n +1)^grp.(1 mod n +1),grp.(2 mod n +1)^grp.(3 mod n +1)]);
	Add(pairs,pairs_group);
od;
pairs_sc := [[AsGroupFRElement(AddingElement(2)),AsGroupFRElement(AddingElement(2))^-1],[AsGroupFRElement(AddingElement(3)),AsGroupFRElement(AddingElement(3))^-1],[AsGroupFRElement(AddingElement(2)),(AsGroupFRElement(AddingElement(2)))^FRElement([[[1],[AddingElement(2)]]],[()],[1])]];


