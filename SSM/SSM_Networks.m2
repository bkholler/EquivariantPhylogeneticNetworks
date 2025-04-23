load("../PhylogeneticNetworks.m2")


-- generate all leaf-labellings for the SSM model on n leaves
ssmStates = n -> (

  states := toList(toList(n:{0,0})..toList(n:{1,1}));

  return delete(null, apply(states, state -> if sum(apply(state, ind -> ind_0)) % 2 == 0 then state))
  )


-- standard ring of probability coordinates
pRing = method(Options => {})
pRing (Number, Number) := Ring => opts -> (k, n) -> QQ[p_(n:0)..p_(n:k-1)]


ssmParameterRing = method(Options => {rVariableName => "r", mVariableName => "m"})
ssmParameterRing Digraph := Ring => opts -> N -> (

	indSet := flatten for e in edges(N) list flatten flatten for i from 0 to 1 list for j from 0 to 1 list for l from 0 to 1 list (toList(e), i, j, l);

	r := getSymbol opts.rVariableName; 
	m := getSymbol opts.mVariableName;

	return R := QQ[s, toList(r_0..r_1) | apply(indSet, i -> m_i)];
	);


ssmReParameterRing = method(Options => {rVariableName => "r", mVariableName => "m", CoefficientRing => QQ})
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

	KK := opts.CoefficientRing;

	return R := KK[s, toList(r_0..r_1) | apply(indSet|condDist, i -> m_i)];
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
ssmNetworkReParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true, CoefficientRing => QQ})
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
	KK := opts.CoefficientRing;
	R := ssmReParameterRing(reticulationEdges, N, CoefficientRing => KK);
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


ssmFourierParameterRing = method(Options => {aVariableName => "a", CoefficientRing => QQ})
ssmFourierParameterRing (List, Digraph) := Ring => opts -> (reticulationEdges, N) -> (

  	int := sort internalVertices(N);
	L := sort delete(null, vertices(N) / (i -> if degree(N, i) == 1 then i));
	n := #L;
	retVert := (reticulationEdges / set // intersect // toList )_0;
	retLeaf := (toList children(N, retVert))_0;
	retParents := sort toList parents(N, retVert);
	T := deleteEdges(N, reticulationEdges|{{retVert, retLeaf}});
	states := {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
	int = delete(retVert, int);

  	a := getSymbol opts.aVariableName;
  	KK := opts.CoefficientRing;

  	return KK[s, toList(r_0..r_1), flatten for e in edges(T) list a_(e, 0, 0, 0)..a_(e, 1, 1, 1), a_({retParents_0, retLeaf}, 0, 0, 0)..a_({retParents_0, retLeaf}, 1, 1, 1), a_({retParents_1, retLeaf}, 0, 0, 0)..a_({retParents_1, retLeaf}, 1, 1, 1)];
)


ssmFourierNetworkParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true, CoefficientRing => QQ})
ssmFourierNetworkParametrization (List, Digraph) := RingMap => opts -> (reticulationEdges, N) -> (

	int := sort internalVertices(N);
	L := sort delete(null, vertices(N) / (i -> if degree(N, i) == 1 then i));
	n := #L;
	retVert := (reticulationEdges / set // intersect // toList )_0;
	retLeaf := (toList children(N, retVert))_0;
	retParents := sort toList parents(N, retVert);
	(u, v) = toSequence retParents;
	rho := (delete(null, apply(int, i -> if degreeIn(N, i) == 0 then i)))_0;
	T := deleteEdges(N, reticulationEdges|{{retVert, retLeaf}});

	-- compute the standard parametrization
	KK := opts.CoefficientRing;
	phi := ssmNetworkReParametrization(reticulationEdges, N, UseStochasticParameters => false, CoefficientRing => KK);

	-- compute the change of coordinates on parameter space
	states := {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
	probToFourierIms := new MutableHashTable;
	for e in edges(T)|{{u, retLeaf}, {v, retLeaf}} do(
		for l1 from 0 to 3 do for l2 from 0 to 3 do(

			(j1, i1) = toSequence states_l1;
			(j2, i2) = toSequence states_l2;

			probToFourierIms#(e, j1, j2, i1, i2) = 1/2*sum(flatten for k1 in {0, 1} list for k2 in {0, 1} list ((-1)^(k1*j1 + k2*j2))*(m_(e, (k1+k2) % 2, i1, i2)))
		);
	);

	R := ssmFourierParameterRing(reticulationEdges, N, CoefficientRing => KK);
	R' := target phi;

	transMatrixIms := drop(gens R, 3) / (i -> last baseName(i)) / (i -> probToFourierIms#(i_0, i_1, i_1, i_2, i_3));
	fourierToTransition := map(R', R, apply({s, r_0, r_1}, i -> sub(i, R')) | transMatrixIms);
	paramL := jacobian matrix fourierToTransition;
	invL := inverse(paramL);
	transitionToFourier := map(R, R', transpose(sub(invL, R)*(transpose matrix {gens R})));

	-- compute the fourier transform of the probabilities
	S' := QQ[p_(toList(n:{0,0}))..p_(toList(n:{1,1}))];

	probL := for x in gens(S') list(

		ind := last baseName x;

		j := apply(ind, i -> i_0);
		i := apply(ind, i -> i_1);

		sum for k in toList((n:0)..(n:1)) list(

			sign := (-1)^(sum toList apply(j, k, (jl, kl) -> jl*kl));
			pInd := apply(k, i, (kl, il) -> {kl, il});

			sign*p_pInd
		)
	);


	(M, H) := (coefficients matrix {probL});
	H = transpose sub(H, KK);
	fourierImages := delete(0_(R'), flatten entries(H*transpose(matrix phi))) / transitionToFourier;
	S := if opts.SourceRing === null then KK[ssmStates(n) / (i -> q_i)] else opts.SourceRing;

	if not opts.UseStochasticParameters then return map(R, S, fourierImages);

	subRules := {r_0 => 1}|(flatten for e in edges(T)|{{u, retLeaf}, {v, retLeaf}} list {a_(e, 0, 0, 1) => 1 - a_(e, 0, 0, 0), a_(e, 0, 1, 1) => 1 - a_(e, 0, 1, 0)});
  	fourierImages = apply(fourierImages, i -> s*sub(i, subRules));
	
	return map(R, S, fourierImages);
)