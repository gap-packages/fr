if fail = LoadPackage("AutoDoc", ">= 2022.07.10") then
    Error("AutoDoc 2022.07.10 or newer is required");
fi;
AutoDoc(rec(
    gapdoc := rec(
        main:="fr.xml",
        files:=["PackageInfo.g"],
    )
));

QUIT;
