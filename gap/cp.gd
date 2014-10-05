DeclareOperation("IteratedOrbit", [ IsFRElement, IsObject]);
DeclareAttribute("OrbitSignalizer", IsFRElement);
################### Branch Groups #############################
DeclareGlobalVariable("START_CP_BRANCH@");
DeclareGlobalFunction("GrigorchukConjugateBranchInit");
DeclareOperation("InitConjugateForBranchGroups", [IsFRGroup,IsList,IsList]);

