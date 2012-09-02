LoadPackage("fr");
SetInfoLevel(InfoFR,1);
dirs := DirectoriesPackageLibrary("fr","tst");
ReadTest(Filename(dirs,"chapter-12.tst"));
ReadTest(Filename(dirs,"chapter-3.tst"));
ReadTest(Filename(dirs,"chapter-4.tst"));
ReadTest(Filename(dirs,"chapter-5-a.tst"));
ReadTest(Filename(dirs,"chapter-5-b.tst"));
ReadTest(Filename(dirs,"chapter-9-a.tst"));
if @FR.dll then
    ReadTest(Filename(dirs,"chapter-9-b.tst"));
else
    Info(InfoFR,1,"The DLL was not compiled, I'm skipping chapter-9-b.tst");
fi;
# hurwitzdirs := DirectoriesPackageLibrary("fr","hurwitz/tst");
# ReadTest( Filename( hurwitzdirs,"utils.tst") );
# ReadTest( Filename( hurwitzdirs,"padicLift.tst") );
# ReadTest( Filename( hurwitzdirs,"hurwitz.tst") );

