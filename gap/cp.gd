#############################################################################
##
#W cp.gd                                                       Thorsten Groth
##
#Y Copyright (C) 2014, Thorsten Groth
##
#############################################################################
##
##  This file declares some attributes for the conjugacy problem
##
#############################################################################

## <#GAPDoc Label="CP">
## <ManSection>
##   <Attr Name="OrbitSignalizer" Arg="g"/>
##   <Returns>The Orbit Signalizer of the group element <A>g</A></Returns>
##   <Description>
##   This attribute computes the orbit signalizer of an element. This is the set
##	 <M>OS(g) := \{g^{|Orb_g(v)|}@v \mid v \in X^*\}</M> where <M>X</M> is the
##	 alphabet of the element <A>g</A> and <M>Orb_g(v)</M> is the orbit of <M>v</M> 
##	 under <M>\langle g \rangle</M>.
## <Example><![CDATA[
## gap> a := MealyElement([[2,2],[2,2]],[(1,2),()],1);
## <Mealy element on alphabet [ 1 .. 2 ] with 2 states>
## gap> OrbitSignalizer(a);
## [ <Mealy element on alphabet [ 1 .. 2 ] with 2 states>, <Trivial Mealy element on alphabet [ 1 .. 2 ]> ]
## ]]></Example>
##   </Description>
## </ManSection>
DeclareAttribute("OrbitSignalizer", IsFRElement);
##
## <ManSection>
##   <Attr Name="FRConjugacyAlgorithm" Arg="G"/>
##   <Returns>A function which solves the conjugacy problem for <A>G</A></Returns>
##   <Description>
##   This attribute stores a function in three arguments which computes a representative
##	 conjugator if exists or fail otherwise. 
##
##   <P/>This attribute is not meant to have a standard setter but to be set if a specialized 
##	 conjugacy algorithm for a certain group is discovered.
## <Example><![CDATA[
## gap> f := FRConjugacyAlgorithm(GrigorchukGroup);
## function( G, g, h ) ... end
## gap> AssignGeneratorVariables(GrigorchukGroup);
## #I  Assigned the global variables [ "a", "b", "c", "d" ]
## gap> f(GrigorchukGroup,a,a^b);
## <Mealy element on alphabet [ 1 .. 2 ] with 5 states>
## ]]></Example>
##   </Description>
## </ManSection>
DeclareAttribute("FRConjugacyAlgorithm", IsFRGroup,2);
##
## <ManSection>
##   <Attr Name="FRBranchGroupConjugacyData" Arg="G"/>
##   <Returns>The initial data for the branch algorithm for <A>G</A></Returns>
##   <Description>
##   This attribute records the data for the branch algorithm. The record has the following components:
##	 <List>
##	  <Mark>initial_conj_dic:</Mark><Item>Dictionary of already known conjugacy pairs with corresponding conjugator tuples.
##			This has to cover at least the TorsionNucleus of <A>G</A></Item>
##	  <Mark>Branchstructure</Mark><Item>Usally calculated by the function BranchStructure</Item>
##	  <Mark>RepSystem</Mark><Item>List of representatives of <M>G/K</M> where <M>K</M> is the branching subgroup of <A>G</A></Item>
##	 </List>
##   </Description>
## </ManSection>
DeclareAttribute("FRBranchGroupConjugacyData", IsFRGroup,"mutable");
## <#/GAPDoc>
