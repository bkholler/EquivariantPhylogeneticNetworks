load("../PhylogeneticNetworks.m2")

pRing = method(Options => {})
pRing (Number, Number) := Ring => opts -> (k, n) -> QQ[p_(n:0)..p_(n:k-1)]


ssmParameterRing = method(Options => {rVariableName => "r", mVariableName => "m"})
ssmParameterRing Digraph := Ring => opts -> N -> (

	indSet := flatten for e in edges(N) list flatten flatten for i from 0 to 1 list for j from 0 to 1 list for l from 0 to 1 list (toList(e), i, j, l);

	r := getSymbol opts.rVariableName; 
	m := getSymbol opts.mVariableName;

	return R := QQ[s, toList(r_0..r_1) | apply(indSet, i -> m_i)];
	);


ssmReParameterRing = method(Options => {rVariableName => "r", mVariableName => "m"})
ssmReParameterRing (List, Digraph) := Ring => opts -> (reticulationEdges, N) -> (

	int := sort internalVertices(N);
	L := sort delete(null, vertices(N) / (i -> if degree(N, i) == 1 then i));
	n := #L;
	retVert := (reticulationEdges / set // intersect // toList )_0;
	retLeaf := (toList children(N, retVert))_0;
	retParents := sort toList parents(N, retVert);
  	T := deleteEdges(N, reticulationEdges|{{retVert, retLeaf}});

	indSet := flatten for e in edges(T) list flatten flatten for i from 0 to 1 list for j from 0 to 1 list for l from 0 to 1 list (toList(e), i, j, l);
	condDist := flatten for e in {{retParents_0, retLeaf}, {retParents_1, retLeaf}} list flatten flatten for i from 0 to 1 list for j from 0 to 1 list for l from 0 to 1 list (toList(e), i, j, l);

	r := getSymbol opts.rVariableName; 
	m := getSymbol opts.mVariableName;

	return R := QQ[s, toList(r_0..r_1) | apply(indSet|condDist, i -> m_i)];
	);


makeSSMHash = (N, R) -> (

	edgeHash := new MutableHashTable;

	for e in edges(T) do(

		M0 := matrix {{m_(e, 0, 0, 0), m_(e, 0, 0, 1)}, {m_(e, 0, 1, 0), m_(e, 0, 1, 1)}};
		M1 := matrix {{m_(e, 1, 0, 0), m_(e, 1, 0, 1)}, {m_(e, 1, 1, 0), m_(e, 1, 1, 1)}};

		edgeHash#e = (M0|M1)||(M1|M0);
	);

	return edgeHash
)


makeSSMReHash = (reticulationEdges, N, R) -> (

	int := sort internalVertices(N);
	L := sort delete(null, vertices(N) / (i -> if degree(N, i) == 1 then i));
	n := #L;
	retVert := (reticulationEdges / set // intersect // toList )_0;
	retLeaf := (toList children(N, retVert))_0;
	retParents := sort toList parents(N, retVert);
  	T := deleteEdges(N, reticulationEdges|{{retVert, retLeaf}});

	edgeHash := new MutableHashTable;

	for e in edges(T) do(

		M0 := matrix {{m_(e, 0, 0, 0), m_(e, 0, 0, 1)}, {m_(e, 0, 1, 0), m_(e, 0, 1, 1)}};
		M1 := matrix {{m_(e, 1, 0, 0), m_(e, 1, 0, 1)}, {m_(e, 1, 1, 0), m_(e, 1, 1, 1)}};

		edgeHash#e = (M0|M1)||(M1|M0);
	);


	for e in apply(retParents, i -> {i, retLeaf}) do(

		M0 := matrix {{m_(e, 0, 0, 0), m_(e, 0, 0, 1)}, {m_(e, 0, 1, 0), m_(e, 0, 1, 1)}};
		M1 := matrix {{m_(e, 1, 0, 0), m_(e, 1, 0, 1)}, {m_(e, 1, 1, 0), m_(e, 1, 1, 1)}};

		edgeHash#e = (M0|M1)||(M1|M0);
	);

	return edgeHash
)


-- k, the number of states for the general Markov model
-- rho, a vertex of T, representing the root of T
-- T, a digraph, representing a phylogenetic tree
ssmTreeParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true})
ssmTreeParametrization Digraph := RingMap => opts -> T -> (

	-- compute internal vertices, leaves, and root
	int := sort internalVertices(T);
	L := sort leaves(graph(graph(T)));
	rho := (delete(null, apply(int, i -> if degreeIn(T, i) == 0 then i)))_0;
	n := #L;
	

	-- make source and target rings
	R := ssmParameterRing(T);
	M := makeSSMHash(T, R);
	k := numcols (values M)_0;
	S := if opts.SourceRing === null then pRing(k, n) else opts.SourceRing;

	-- for each state at the leaves, compute the probability by marginalizing out over all internal states
	phi := for leafState in toList((n:0)..(n:k-1)) list(

			sum for intState in toList(toList(#int:0)..toList(#int:k-1)) list(

				states := hashTable(apply(L, i -> i => leafState_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));

				(r_((states#rho) % 2))*product(apply(edges(T), e -> (M#e)_(states#(e_0), states#(e_1))))
				)
			);


	if not opts.UseStochasticParameters then return map(R, S, phi);

	subRules := flatten for e in edges(T) list(

		for i from 0 to 1 list(

			(M#e)_(i, k-1) => (1 - sum apply(k-1, j -> (M#e)_(i, j)))
		)
	);

	subRules = subRules|{r_1 => 1/2 - r_0};

	phi = apply(phi, i -> s*sub(i, subRules));
	
	return map(R, S, phi)
	);


-- k, the number of states for the general Markov model
-- N, a digraph, representing a phylogenetic network
-- reticulationEdges, a list of edges of N, in the form {i, j}, representing the reticulation edges of N
ssmNetworkParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true})
ssmNetworkParametrization (List, Digraph) := RingMap => opts -> (reticulationEdges, N) -> (

	-- compute internal vertices, leaves, and root
	-- also compute underlying trees
	T1 := deleteEdges(N, {reticulationEdges_0});
	T2 := deleteEdges(N, {reticulationEdges_1});
	int := sort internalVertices(T1);
	L := sort leaves(graph(graph(T1)));
	rho := (delete(null, apply(int, i -> if degreeIn(N, i) == 0 then i)))_0;
	n := #L;

	-- make source and target rings
	R := ssmParameterRing(N);
	M := makeSSMHash(N, R);
	k := numcols (values M)_0;
	S := if opts.SourceRing === null then pRing(k, n) else opts.SourceRing;

	phi := for leafState in toList((n:0)..(n:k-1)) list(

			sum for intState in toList((#int:0)..(#int:k-1)) list(

				states := hashTable(apply(L, leafState, (i,j) -> i => j) | apply(int, intState, (i,j) -> i =>j));

				(r_((states#rho) % 2))*(t*product(apply(edges(T1), e -> (M#e)_(states#(e_0), states#(e_1)))) + (1-t)*product(apply(edges(T2), e -> (M#e)_(states#(e_0), states#(e_1)))))
				)
			);

	if not opts.UseStochasticParameters then return map(R, S, phi);

	subRules := flatten for e in edges(N) list(

		for i from 0 to 1 list(

			(M#e)_(i, k-1) => (1 - sum apply(k-1, j -> (M#e)_(i, j)))
		)
	);

	subRules = subRules|{r_1 => 1/2 - r_0};

	phi = apply(phi, i -> s*sub(i, subRules));
	
	return map(R, S, phi)
	);


-- k, the number of states for the general Markov model
-- N, a digraph, representing a phylogenetic network
-- reticulationEdges, a list of edges of N, in the form {i, j}, representing the reticulation edges of N
ssmNetworkReParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true})
ssmNetworkReParametrization (List, Digraph) := RingMap => opts -> (reticulationEdges, N) -> (

	-- compute internal vertices, leaves, and root
	-- also compute underlying trees
	int := sort internalVertices(N);
	L := sort delete(null, vertices(N) / (i -> if degree(N, i) == 1 then i));
	n := #L;
	retVert := (reticulationEdges / set // intersect // toList )_0;
	retLeaf := (toList children(N, retVert))_0;
	retParents := sort toList parents(N, retVert);
	(u, v) = toSequence retParents;
	rho := (delete(null, apply(int, i -> if degreeIn(N, i) == 0 then i)))_0;
  	T := deleteEdges(N, reticulationEdges|{{retVert, retLeaf}});
	int = delete(retVert, int);

	-- make source and target rings
	R := ssmReParameterRing(reticulationEdges, N);
	M := makeSSMReHash(reticulationEdges, N, R);
	k := numcols (values M)_0;
	S := if opts.SourceRing === null then pRing(k, n) else opts.SourceRing;

	phi := for leafState in toList((n:0)..(n:k-1)) list(

			sum for intState in toList((#int:0)..(#int:k-1)) list(

				states := hashTable(apply(L, leafState, (i,j) -> i => j) | apply(int, intState, (i,j) -> i =>j));

				
				(r_((states#rho) % 2))*product(apply(edges(T), e -> (M#e)_(states#(e_0), states#(e_1))))*((M#{u, retLeaf})_(states#u, states#retLeaf) + (M#{v, retLeaf})_(states#v, states#retLeaf))
				)
			);

	if not opts.UseStochasticParameters then return map(R, S, phi);

	subRules := flatten for e in edges(T)|apply(retParents, i -> {i, retLeaf}) list(

		for i from 0 to 1 list(

			(M#e)_(i, 3) => (1 - sum apply(3, j -> (M#e)_(i, j)))
		)
	);

	subRules = subRules|{r_1 => 1/2 - r_0};

	phi = apply(phi, i -> s*sub(i, subRules));
	
	return map(R, S, phi)
	);
