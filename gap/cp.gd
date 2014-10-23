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
##   <Returns>...</Returns>
##   <Description>
##     ...
## <Example><![CDATA[
## gap> OrbitSignalizer(...);
## ...
## ]]></Example>
##   </Description>
## </ManSection>
DeclareAttribute("OrbitSignalizer", IsFRElement);
##
## <ManSection>
### ...
## </ManSection>
DeclareAttribute("FRConjugacyAlgorithm", IsFRGroup,2);
##
## <ManSection>
### ...
## </ManSection>
DeclareAttribute("FRBranchGroupConjugacyData", IsFRGroup,"mutable");
## <#/GAPDoc>
