load("../PhylogeneticNetworks.m2")
load("../MatroidDistinguish.m2")
-- This loads our new package which contains our implementation of Algorithm 1
needsPackage "MultigradedImplicitization"

-- This loads the PhylogeneticTrees package for macaulay2
needsPackage "PhylogeneticTrees"
needsPackage "Permutations"


-- n, the number of leaves
-- M, a model, see PhylogeneticTrees.m2
-- Outputs the parameterization of a cyclically labeled n-leaf sunlet network with the leaf labeled 1 as the reticulation vertex
-- The output is given as a list whose entries correspond to the fourier coordinates in lexicographic order on the label sequence
sunletParam = (n, M) -> (

	-- This makes 
	indS := leafColorings(n,M);
	GH := new MutableHashTable from apply(group(M), toList(0..#group(M)-1), (i, j) -> i => j);

	-- This makes the ring of parameters which is the codomain of phi
	-- the "a" parameters correspond to the leaves of the network while the "b" parameters correspond to the internal edges
	indR := flatten for i from 1 to n list apply(group(M), j -> {i, GH#j});
	R := QQ[apply(indR, k -> a_k) | apply(indR, k -> b_k)];
	

	images := for g in indS list(

		aProd := product(apply(n, i -> a_{i+1, GH#(groupElToRep(g_i, M))}));
		bProd1 := product(for j from 0 to #g - 2 list b_{j+1, GH#(groupElToRep(sum(g_(toList(0..j))), M))});
    	bProd2 := product(for j from 1 to #g- 1 list b_{j+1, GH#(groupElToRep(sum(g_(toList(1..j))), M))});
    	aProd*(bProd1 + bProd2)
		);

	return images
	)


oneCycleNetworkParam = (reticulationEdges, N, M) -> (

	indS := leafColorings(n,M);
	GH := new MutableHashTable from apply(group(M), toList(0..#group(M)-1), (i, j) -> i => j);

	T1 := deleteEdges(N, {reticulationEdges_0});
	T2 := deleteEdges(N, {reticulationEdges_1});
	L := sort leaves underlyingGraph T1;
	n := #L;

	-- This makes the ring of parameters which is the codomain of phi
	-- the "a" parameters correspond to the leaves of the network while the "b" parameters correspond to the internal edges
	indR = flatten for i from 0 to #edges(N)-1 list apply(group(M), j -> {i, GH#j});
	R := QQ[apply(indR, k -> a_k)];

	H1 := hashTable apply(edges T1, e -> {e, select(first connectedComponents underlyingGraph deleteEdges(T1, {e}), l -> member(l, L))});
	H2 := hashTable apply(edges T2, e -> {e, select(first connectedComponents underlyingGraph deleteEdges(T2, {e}), l -> member(l, L))});
	H1 = applyPairs(H1, (k, v) -> (position(edges(N), j -> j == k), v / (i -> i - 1)));
	H2 = applyPairs(H2, (k, v) -> (position(edges(N), j -> j == k), v / (i -> i - 1)));

	images := for g in indS list(

		term1 := product apply(keys H1, i -> a_{i, GH#(groupElToRep((sum(g_(H1#i))), M))});
		term2 := product apply(keys H2, i -> a_{i, GH#(groupElToRep((sum(g_(H2#i))), M))});
		term1 + term2
	);

	return images
)

groupElToRep = (g, M) -> (

	B := buckets(M);
	return first flatten select(B, b -> member(g, b))
)


end