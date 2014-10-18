DeclareOperation("IteratedOrbit", [ IsFRElement, IsObject]);
DeclareAttribute("OrbitSignalizer", IsFRElement);
################### Groups with special algorithm   ###########
DeclareAttribute("FRConjugacyAlgorithm", IsFRGroup);
################### Branch Groups #############################
DeclareGlobalVariable("START_CP_BRANCH@");
DeclareGlobalFunction("GrigorchukConjugateBranchInit");
DeclareGlobalFunction("GuptaSidkiConjugateBranchInit");
DeclareOperation("InitConjugateForBranchGroups", [IsFRGroup,IsList,IsList]);

