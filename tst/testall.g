LoadPackage("FR");
SetInfoLevel(InfoFR,1);
dirs := DirectoriesPackageLibrary("FR","tst");
TestDirectory(dirs, rec(exitGAP := true));
FORCE_QUIT_GAP(1);
