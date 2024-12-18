load("../PhylogeneticNetworks.m2")


-- N, a phylogenetic network
-- Creates the ring of parameters whose variables are the root parameters and entries of the transition matrices
gmmParameterRing = method(Options => {rVariableName => "r", mVariableName => "m"})
gmmParameterRing (Number, Digraph) := Ring => opts -> (k, N) -> (

	indSet := flatten for e in edges(N) list flatten for i from 0 to k-1 list for j from 0 to k-1 list (toList(e), i, j);

	r := getSymbol opts.rVariableName; 
	m := getSymbol opts.mVariableName;

	return R := QQ[apply(k, i -> r_i) | apply(indSet, i -> m_i)];
	);


-- k, number of states of the GMM model
-- n, number of leaves of the phylogenetic network on which the model is to be defined
pRing = method()
pRing (Number, Number) := Ring => opts -> (k, n) -> QQ[p_(n:0)..p_(n:k-1)]


-- k, the number of states for the general Markov model
-- rho, a vertex of T, representing the root of T
-- T, a digraph, representing a phylogenetic tree
gmmTreeParametrization = method(Options => {SourceRing => null})
gmmTreeParametrization (Number, Digraph) := RingMap => opts -> (k, T) -> (

	-- compute internal vertices, leaves, and root
	int := sort internalVertices(T);
	L := sort leaves(graph(graph(T)));
	rho := (delete(null, apply(int, i -> if degreeIn(T, i) == 0 then i)))_0;
	n := #L;

	-- make source and target rings
	R := gmmParameterRing(k, T);
	S := if opts.SourceRing === null then pRing(k, n) else opts.SourceRing;

	-- for each state at the leaves, compute the probability by marginalizing out over all internal states
	phi := for leafState in toList((n:0)..(n:k-1)) list(

			sum for intState in toList(toList(#int:0)..toList(#int:k-1)) list(

				states := hashTable(apply(L, i -> i => leafState_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));

				(r_(states#rho))*product(apply(edges(T), e -> (m_(e, states#(e_0), states#(e_1)))))
				)
			);

	return map(R, S, phi)
	);


-- k, the number of states for the general Markov model
-- N, a digraph, representing a phylogenetic network
-- reticulationEdges, a list of edges of N, in the form {i, j}, representing the reticulation edges of N
gmmNetworkParametrization = method(Options => {SourceRing => null})
gmmNetworkParametrization (Number, List, Digraph) := RingMap => opts -> (k, reticulationEdges, N) -> (

	-- compute internal vertices, leaves, and root
	-- also compute underlying trees
	T1 := deleteEdges(N, {reticulationEdges_0});
	T2 := deleteEdges(N, {reticulationEdges_1});
	int := sort internalVertices(T1);
	L := sort leaves(graph(graph(T1)));
	root := (delete(null, apply(int, i -> if degreeIn(N, i) == 0 then i)))_0;
	n := #L;

	-- make source and target rings
	R := gmmParameterRing(k, N);
	S := if opts.SourceRing === null then pRing(k, n) else opts.SourceRing;

	phi := for leafState in toList((n:0)..(n:k-1)) list(

			sum for intState in toList((#int:0)..(#int:k-1)) list(

				states := hashTable(apply(L, leafState, (i,j) -> i => j) | apply(int, intState, (i,j) -> i =>j));

				(r_(states#rho))*(product(apply(edges(T1), e -> m_(e, states#(e_0), states#(e_1))  )) + product(apply(edges(T2), e -> m_(e, states#(e_0), states#(e_1)) )))
				)
			);

	return map(R, S, phi)
	);


-- A,B, a pair of list which form a bipartition of [n] representing a split
-- S, the ring of the phylogenetic Markov model on a n-leaf network with k states
-- returns the flattening matrix corresponding to the partition A|B
flat = method(Options => {})
flat (List, List, Ring) := Matrix => opts -> (A, B, S) -> (
	
	-- compute k and n from the ring S
	k := (last baseName (last gens S))_0;
	n := #(last baseName S_0);

	-- make the set of all joint states of the leaves in A and leaves in B
	indA := toList(#A:0)..toList(#A:(k- 1));
	indB := toList(#B:0)..toList(#B:(k- 1));
	x := first baseName S_0;

	-- build the flattening matring
	return matrix for a in indA list(

		for b in indB list(

			curInd := new MutableHashTable from apply(#A, i -> A_i => a_i) | apply(#B, i -> B_i => b_i);

			(x_(toSequence apply(n, i -> curInd#(i+1))))_S
			)
		);
	);
