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

DeclareCategory("IsFRBiset",
        IsAssociativeROpDProd and IsAssociativeLOpDProd);

DeclareRepresentation("IsFRBisetByFRMachine", IsFRBiset, []);

BindGlobal("FRBISET_FAMILY",
        NewFamily("FRBisetsFamily",IsFRBiset));

DeclareOperation("BisetByFRMachine", [IsFRMachine]);
DeclareOperation("BisetByFRSemigroup", [IsFRSemigroup]);
DeclareOperation("BisetByHomomorphism", [IsMagmaHomomorphism]);

DeclareProperty("IsLeftFree", IsFRBiset);
DeclareProperty("IsRightFree", IsFRBiset);
DeclareProperty("IsLeftTransitive", IsFRBiset);
DeclareProperty("IsRightTransitive", IsFRBiset);

#E bisets.gd. . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
