needsPackage "Graphs";

-- T, a tree in the form of a directed graph
-- returns all internal vertices of T
internalVertices = T -> delete(null, vertices(T) / (i -> if degree(T, i) > 1 then i))


-- T, a tree in the form of a directed graph
-- L, a subset of the leaves of T
-- i, an internal node
-- returns all leaves in L which are descendants of i
leafDescendants = (i, L, T) -> (

  desc := descendants(T, i);

  return delete(null, apply(L, i -> if member(i, desc) then i))
  )


-- generates all 12 distinct 4-cycle networks
generate4Cycles = () -> (

	N := digraph(toList(1..8), {{5,1}, {6,2}, {7,3}, {8,4}, {6,5}, {8,5}, {7,8}, {7,6}});

	adj := adjacencyMatrix N;

	perms := permutations({0,1,2,3});

	permsModLeafSwitch := while #perms > 0 list(

		sigma := perms_0;

		perms = delete(sigma_{0,3,2,1}, perms);
		perms = delete(sigma, perms);

		
		sigma
		);

	fourLeafNets := for sigma in permsModLeafSwitch list digraph(toList(1..8), adj_(sigma|{4,5,6,7}));

	return fourLeafNets
	)

-- generates all 12 distinct 4-cycle networks and returns a pair consisting of the network and its corresponding permutation
generate4CyclesWithPerms = () -> (

	N := digraph(toList(1..8), {{5,1}, {6,2}, {7,3}, {8,4}, {6,5}, {8,5}, {7,8}, {7,6}});

	adj := adjacencyMatrix N;

	perms := permutations({0,1,2,3});

	permsModLeafSwitch := while #perms > 0 list(

		sigma := perms_0;

		perms = delete(sigma_{0,3,2,1}, perms);
		perms = delete(sigma, perms);

		
		sigma
		);

	fourLeafNets := for sigma in permsModLeafSwitch list digraph(toList(1..8), adj_(sigma|{4,5,6,7}));

	return fourLeafNets
	)


-- t, an integer
-- f, a polynomial in the source of phi
-- phi, a ring map 
-- samples t points from the image of phi over a finite field and evaluates f on all points. 
-- if f ever evaluates to a nonzero number then returns false else it returns true
probabilisticMembershipTest = (t, f, phi) -> (


	FF := ZZ/nextPrime(32003);

	for t from 1 to t do(

		params := apply(gens target phi, i -> i => random(FF));
		pt := apply(gens source phi, j -> j => sub(phi(j), params));

		if not sub(f, pt) == 0_FF then return false;
		);
	
	return true;
	)


-- builds a polynomial in the minors of flat1 and flat2 correspond to inds with coefficients given by coeffs
fillTemplate = (inds, coeffs, flat1, flat2) -> (

  return sum for i from 0 to #inds - 1 list(
	
    ind  = inds_i;
    (R1, C1) = ind_0;
    (R2, C2)= ind_1;

    (coeffs_i)*(det(flat1_C1^R1))*(det(flat2_C2^R2))
    );
  ); 