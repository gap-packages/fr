#############################################################################
##
#W init.g                                                   Laurent Bartholdi
##
#Y Copyright (C) 2006-2016, Laurent Bartholdi
##
#############################################################################
##
##  This file reads the declarations of the packages' new objects
##
#############################################################################

#I introducing globally the NC versions of PreImages...
if not IsBound( PreImagesNC ) then
    BindGlobal( "PreImagesNC", PreImages );
fi;
if not IsBound( PreImagesRepresentativeNC ) then
    BindGlobal( "PreImagesRepresentativeNC", PreImagesRepresentative );
fi;

#############################################################################
BindGlobal("Fr", rec());

Fr.POSTHOOK := []; # to be processed at the end

Fr.locals := rec(); # a record to store locals in the package

#############################################################################
##
#I Create info class to be able to debug loading
##
DeclareInfoClass("InfoFR");
SetInfoLevel(InfoFR, 1);
#############################################################################

#############################################################################
##
#R Read the declaration files.
##
ReadPackage("fr", "gap/helpers.gd");
ReadPackage("fr", "gap/perlist.gd");
ReadPackage("fr", "gap/frmachine.gd");
ReadPackage("fr", "gap/frelement.gd");
ReadPackage("fr", "gap/mealy.gd");
ReadPackage("fr", "gap/group.gd");
ReadPackage("fr", "gap/vector.gd");
ReadPackage("fr", "gap/algebra.gd");
ReadPackage("fr", "gap/bisets.gd");
ReadPackage("fr", "gap/examples.gd");
ReadPackage("fr", "gap/cp.gd");

Fr.locals.nql := IsBound(IsLpGroup);

if not Fr.locals.nql then # shut up warnings in case LpGroups is not present
    Perform(["IsLpGroup","IsElementOfLpGroup","LPresentedGroup",
            "ElementOfLpGroup","SetEmbeddingOfAscendingSubgroup"], function(w)
        BindGlobal(w, fail);
        Add(Fr.POSTHOOK,function() MakeReadWriteGlobal(w); UnbindGlobal(w); end);
    end);
fi;
