#if fail = LoadPackage("AutoDoc", ">= 2016.01.21") then
#    Error("AutoDoc 2016.01.21 or newer is required");
#fi;
#AutoDoc(rec(gapdoc := rec(files:=["PackageInfo.g"])));

MakeGAPDocDoc("doc","fr",
  ["../gap/frmachine.gd","../gap/frelement.gd","../gap/mealy.gd",
   "../gap/group.gd","../gap/vector.gd","../gap/algebra.gd",
   "../gap/examples.gd","../gap/helpers.gd","../gap/perlist.gd","../gap/cp.gd",
   "../PackageInfo.g"],"fr","../../..");
CopyHTMLStyleFiles("doc");
GAPDocManualLab("fr");

QUIT;
