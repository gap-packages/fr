LieAlgebraForgetfulFunctor := function(a)
    local l;
    if IsAlgebra(a) and IsAssociative(a) then
        l := AlgebraByGenerators(LeftActingDomain(a),List(GeneratorsOfAlgebra(a),LieObject));
        return AlgebraHomomorphismByFunction(a,l,LieObject);
    elif IsAlgebraHomomorphism(a) then
        l := List([Source(a),Range(a)],LieAlgebraForgetfulFunctor);
        return AlgebraHomomorphismByFunction(Range(l[1]),Range(l[2]),x->LieObject(Image(a,x![1])));
    else
        TryNextMethod();
    fi;
end;