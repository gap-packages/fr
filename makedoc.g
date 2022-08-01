if fail = LoadPackage("AutoDoc", ">= 2019.04.10") then
    Error("AutoDoc 2019.04.10 or newer is required");
fi;
AutoDoc(rec(
    scaffold:=rec(
        LaTeXOptions := rec( ExtraPreamble := """
            \usepackage{graphicx}
        """ ),
        includes:=["fr.xml"],
        bib := "frbib.xml",
    ),
    gapdoc := rec(
        main:="_main.xml",
        files:=["PackageInfo.g"],
    )
));

QUIT;
