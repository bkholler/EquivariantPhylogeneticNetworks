load("PhylogeneticNetworks")


-- N, a digraph, representing a phylogenetic network
-- sym1, a string,  used for root parameters
-- sym2, a string, is used for transition matrix entries
-- Creates the ring of parameters whose variables are the root parameters and entries of the transition matrices
gmmParamRing = method(Options => {rVariableName => "r", mVariableName => "m"})
gmmParamRing = (Number, Digraph) := Ring => opts -> (k, N) -> (

	indSet := flatten for e in edges(N) list flatten for i from 0 to k-1 list for j from 0 to k-1 list (toList(e), i, j);

	r := getSymbol sym1; 
	m := getSymbol sym2;

	return R := QQ[apply(k, i -> r_i) | apply(indSet, i -> m_i)];
	);


-- k, the number of states for the general Markov model
-- rho, a vertex of T, representing the root of T
-- T, a digraph, representing a phylogenetic tree
gmmTreeParam = (k, rho, T) -> (

	int := sort intVertices(T);
	L := sort leaves(graph(graph(T)));
	n := #L;
	S := gmmParamRing(k, T, "r", "m");

	phi := for leafState in toList((n:0)..(n:k-1)) list(

			sum for intState in toList(toList(#int:0)..toList(#int:k-1)) list(

				states := hashTable(apply(L, i -> i => leafState_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));

				(r_(states#rho))*product(apply(edges(T), e -> (m_(e, states#(e_0), states#(e_1)))))
				)
			);

	return phi
	);


-- k, the number of states for the general Markov model
-- rho, a vertex of T, representing the root of T
-- N, a digraph, representing a phylogenetic network
-- retEdges, a list of edges of N, in the form {i, j}, representing the reticulation edges of N
gmmNetParam = (k, rho, retEdges, N) -> (

	T1 := deleteEdges(N, {retEdges_0});
	T2 := deleteEdges(N, {retEdges_1});
	int := sort intVertices(T1);
	L := sort leaves(graph(graph(T1)));
	n := #L;
	root := (delete(null, apply(int, i -> if degreeIn(N, i) == 0 then i)))_0;
	S := gmmParamRing(k, N, "r", "m");
	use S;

	phi := for leafState in toList((n:0)..(n:k-1)) list(

			sum for intState in toList((#int:0)..(#int:k-1)) list(

				states := hashTable(apply(L, leafState, (i,j) -> i => j) | apply(int, intState, (i,j) -> i =>j));

				(r_(states#rho))*(product(apply(edges(T1), e -> m_(e, states#(e_0), states#(e_1))  )) + product(apply(edges(T2), e -> m_(e, states#(e_0), states#(e_1)) )))
				)
			);

	return phi
	);


-- n, the number of leaves of a network
-- k, the number of states for the general Markov model
-- A,B, a bipartition of [n] representing a split
-- S, the ring of the phylogenetic Markov model on a n-leaf network with k states. 
flat = (k, n, A, B, S) -> (

	indA := toList(#A:0)..toList(#A:(k- 1));
	indB := toList(#B:0)..toList(#B:(k- 1));
	x := first baseName S_0;

	return matrix for a in indA list(

		for b in indB list(

			curInd := new MutableHashTable from apply(#A, i -> A_i => a_i) | apply(#B, i -> B_i => b_i);

			(x_(toSequence apply(n, i -> curInd#(i+1))))_S
			)
		);
	);


