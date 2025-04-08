-- evaluate a matrix at a random point over the field opts.CoefficientRing
specialize = method(Options => {CoefficientRing => ZZ/32003})
specialize Matrix := Matrix => opts -> A -> (

    R := ring A;
    KK := opts.CoefficientRing;
    return sub(A, apply(gens R, i -> i => random(KK)));
)


-- implementation of Algorithm 3.6 from https://arxiv.org/pdf/1909.13754
matroidDistinguish = method(Options => {Dim1 => null, Dim2 => null, MaxEntryDegree => null, Tolerance => 1*10^(-10), Verbose => true})
matroidDistinguish (Number, Number, Matrix, Matrix) := (Boolean, List) => opts -> (t, d, J1, J2) -> (

    -- Pick a coefficient field over which we can do amplification by independent trials
    -- Then specialize the jacobian by plugging in random parameter values over the finite field
    alpha := if opts.MaxEntryDegree === null then max(flatten(apply(flatten entries J1, degree))) else opts.MaxEntryDegree;
    pp := nextPrime(d*alpha + 10^6);
    KK := ZZ/pp;
    l := ceiling(log(2, opts.Tolerance)/log(2, alpha/pp));
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
        
        -- check if the rank differs
        if r1 < r2 then(
            if opts.Verbose then print("Found potential certificate. Beginning verification over finite field.");
            if schwartzZippelCertify(pp, l, r1, J1_S) then return {r1, r2, S} else if opts.Verbose then print("Potential certificate was erroneous. Restarting certificate search.");
        );

        if r1 > r2 then(
            if opts.Verbose then print("Found potential certificate. Beginning verification over finite field.");
            if schwartzZippelCertify(pp, l, r2, J2_S) then return {r1, r2, S} else if opts.Verbose then print("Potential certificate was erroneous. Restarting certificate search.");
        );
    );
)


-- checks if the matrix J has rank at most r, l times over the finite field ZZ.p
schwartzZippelCertify = (p, l, r, J) -> (

    KK := ZZ/p;

    for i from 1 to l do if rank(specialize(J, CoefficientRing => KK)) > r then return false;

    return true
)


-- creates a new matrix which encodes the sparsity pattern of the input matrix
supportMatrix = A -> (

    return matrix for i from 0 to numrows(A)-1 list for j from 0 to numcols(A)-1 list if A_(i, j) != 0 then 1 else 0
)



-- returns the indices of the circuit corresponding to the basis of the matroid of A given by the first r entries of sigma
-- and the additional column given by ind
fundamentalCircuit = (r, sigma, ind, A) -> (

	A' := A_sigma;
	A' = reducedRowEchelonForm(A');
	if rank(submatrix(A', , toList(0..r-1))) != r then error "matrix is not in appropriate form";
	C := submatrix(A', ,toList(0..r-1))|submatrix(A', , {ind});
	K := flatten entries gens ker C;

	return delete(null, apply(drop(K, -1), take(sigma, r), (i, j) -> if i != 0 then j) | {sigma_ind})
)


-- removes all zero rows from the matrix A
nonZeroRows = A -> (

	keepRows := for i from 0 to numrows(A)-1 list if A^{i} != 0 then i else continue;
	return A^keepRows;
)