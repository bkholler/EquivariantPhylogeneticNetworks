specialize = method(Options => {CoefficientRing => ZZ/32003})
specialize Matrix := Matrix => opts -> A -> (

    R := ring A;
    KK := opts.CoefficientRing;
    return sub(A, apply(gens R, i -> i => random(KK)));
)

matroidDistinguish = method(Options => {Dim1 => null, Dim2 => null, MaxEntryDegree => null, Tolerance => 1*10^(-10), Verbose => true})
matroidDistinguish (Number, Number, Matrix, Matrix) := (Boolean, List) => opts -> (t, d, J1, J2) -> (

    -- Pick a coefficient field over which we can do amplification by independent trials
    -- Then specialize the jacobian by plugging in random parameter values over the finite field
    alpha := if opts.MaxEntryDegree === null then max(flatten(apply(flatten entries J1, degree))) else opts.MaxEntryDegree;
    p := nextPrime(d*alpha + 1);
    KK := ZZ/p;
    l := ceiling(log(2, opts.Tolerance)/log(2, alpha/p));
    numJ1 := specialize(J1, CoefficientRing => KK);
    numJ2 := specialize(J2, CoefficientRing => KK);

    -- compute dimensions and throw an error if they are not the same
    dim1 := if opts.Dim1 === null then max(apply(100, i -> rank specialize(J1))) else opts.Dim1;
    dim2 := if opts.Dim2 === null then max(apply(100, i -> rank specialize(J2))) else opts.Dim2;

    if dim1 != dim2 then error("This method is currently only implemented for models of the same dimension");    

    for i from 1 to t do(

        S := take(random(toList(0..numcols(J1)-1)), d);
        r1 := rank(numJ1_S);
        r2 := rank(numJ2_S);
        
        -- check if the rank differes
        if r1 < r2 then(
            if opts.Verbose then print("Found potential certificate. Beginning verification over finite field.");
            if schwartzZippelCertify(p, l, r2, J2_S) then return {r1, r2, S} else if opts.Verbose then print("Potential certificate was erroneous. Restarting certificate search.");
        );

        if r1 < r2 then(
            if opts.Verbose then print("Found potential certificate. Beginning verification over finite field.");
            if schwartzZippelCertify(p, l, r1, J1_S) then return {r1, r2, S} else if opts.Verbose then print("Potential certificate was erroneous. Restarting certificate search.");
        );
    );
)

-- 
schwartzZippelCertify = (p, l, r, J) -> (

    KK := ZZ/p;

    for i from 1 to l do if rank(specialize(J, CoefficientRing => KK)) > r then return false;

    return true
)