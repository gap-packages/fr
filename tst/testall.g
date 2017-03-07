LoadPackage("FR");
SetInfoLevel(InfoFR,1);
dirs := DirectoriesPackageLibrary("FR","tst");
Test(Filename(dirs,"chapter-12.tst"));
Test(Filename(dirs,"chapter-3.tst"));
Test(Filename(dirs,"chapter-4.tst"));
Test(Filename(dirs,"chapter-5-a.tst"));
Test(Filename(dirs,"chapter-5-b.tst"));
