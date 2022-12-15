#############################################################################
##
#W bisets.gd                                                Laurent Bartholdi
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  This file defines general bisets.
##
#############################################################################

# BISET:
# has left/right acting domains, and elementsfamily.
# two representations:
# (1) as homomorphism; then elements are elements in the range
# (2) as machines (maybe twisted, in that the group entries can be from wherever); then elements are pairs (integer,group element) in the "natural" basis
# BASIS: for left-free, stored as (permutation,sequence of group elements).
# OPERATIONS:
# left/right product of elements by acting domains
# tensor product, contragredient
# for a given basis, extract machine, or wreath recursion
# test, and produce, isomorphism of two bisets, in the form of a basis of the second biset that makes them identical.
# multiply biset by homomorphisms on left/right to change coordinates.
# decompose as amalgamated free product, etc.
# congruence of G-G-bisets, by conjugation.

DeclareCategory("IsBiset", IsObject);

DeclareCategory("IsFRBiset",
        IsBiset and IsAssociativeROpDProd and IsAssociativeLOpDProd);

DeclareRepresentation("IsFRBisetByFRMachineRep", IsFRBiset and IsAttributeStoringRep, []);
DeclareRepresentation("IsFRBisetByFRSemigroupRep", IsFRBiset and IsAttributeStoringRep, []);
DeclareRepresentation("IsFRBisetByHomomorphismRep", IsFRBiset and IsAttributeStoringRep, []);

DeclareCategory("IsBisetElement",
        IsObject);

DeclareRepresentation("IsBisetElementByPair", IsBisetElement and IsComponentObjectRep, []);
DeclareRepresentation("IsBisetElementByElement", IsBisetElement and IsComponentObjectRep, []);

BindGlobal("FRBISET_FAMILY",
        NewFamily("FRBisetsFamily",IsObject));

DeclareOperation("BisetByFRMachine", [IsFRMachine]);
DeclareOperation("BisetByFRSemigroup", [IsFRSemigroup]);
DeclareSynonym("BisetByFRMonoid", BisetByFRSemigroup);
DeclareSynonym("BisetByFRGroup", BisetByFRSemigroup);

DeclareAttribute("DualBiset", IsFRBiset);
DeclareOperation("TensorProductOp", [IsList,IsFRBiset]);

DeclareOperation("BisetElement", [IsFRBiset,IsMultiplicativeElement,IsObject]);
DeclareOperation("BisetElement", [IsFRBiset,IsMultiplicativeElement]);
        
DeclareProperty("IsLeftFree", IsFRBiset);
DeclareProperty("IsRightFree", IsFRBiset);
DeclareProperty("IsLeftTransitive", IsFRBiset);
DeclareProperty("IsRightTransitive", IsFRBiset);

DeclareCategory("IsBisetBasis", IsBasis);
DeclareCategory("IsLeftBisetBasis", IsBisetBasis);
DeclareCategory("IsRightBisetBasis", IsBisetBasis);
DeclareAttribute("Basis", IsFRBiset);
DeclareAttribute("LeftBasis", IsFRBiset);
DeclareOperation("LeftBasis", [IsFRBiset,IsList]);
DeclareOperation("LeftBasis", [IsFRBiset,IsPerm,IsList]);
DeclareOperation("LeftBasis", [IsFRBiset,IsTransformation,IsList]);
DeclareOperation("LeftBasis", [IsFRBiset,IsList,IsList]);
DeclareAttribute("RightBasis", IsFRBiset);
DeclareAttribute("CanonicalBasis", IsFRBiset);

DeclareAttribute("WreathRecursion", IsFRBiset);
DeclareOperation("WreathRecursion", [IsFRBiset,IsLeftBisetBasis]);

DeclareAttribute("FRMachineOfBiset", IsFRBiset);
DeclareOperation("FRMachine", [IsFRBiset]);
DeclareOperation("FRMachine", [IsFRBiset,IsLeftBisetBasis]);

