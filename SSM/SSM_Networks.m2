load "../PhylogeneticNetworks.m2"


-- generate all leaf-labellings for the SSM model on n leaves
ssmStates = n -> (

  states := toList(toList(n:{0,0})..toList(n:{1,1}));

  return delete(null, apply(states, state -> if sum(apply(state, ind -> ind_0)) % 2 == 0 then state))
  )

-- N, a digraph representing a phylogenetic network
-- creates the ring of fourier parameters for the SSM model on N
ssmFourierParameterRing = method(Options => {aVariableName => "a", CoefficientRing => QQ})
ssmFourierParameterRing (List, Digraph) := Ring => opts -> (reticulationEdges, N) -> (

  int := sort internalVertices(N);
	L := sort delete(null, vertices(N) / (i -> if degree(N, i) == 1 then i));
	n := #L;
	retVert := (reticulationEdges / set // intersect // toList )_0;
	retLeaf := (toList children(N, retVert))_0;
	retParents := sort toList parents(N, retVert);
  T := deleteEdges(N, reticulationEdges);
  a := getSymbol opts.aVariableName;
  KK := opts.CoefficientRing;

  return KK[toList(r_0..r_1), flatten for e in edges(T) list a_(e, 0, 0, 0)..a_(e, 1, 1, 1), a_({retParents_0, retLeaf}, 0, 0, 0)..a_({retParents_0, retLeaf}, 1, 1, 1), a_({retParents_1, retLeaf}, 0, 0, 0)..a_({retParents_1, retLeaf}, 1, 1, 1)];
  )

-- T, a digraph representing a phylogenetic tree
-- creates the parametrization of the strand symmetric model on T in the fourier coordinates
-- By default the stochastic restrictions are enforced but the parametrization is re-homogenized 
ssmFourierTreeParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true, CoefficientRing => QQ})
ssmFourierTreeParametrization Digraph := RingMap => opts -> T -> (

  -- compute internal vertices, leaves, and displayed trees
  int := sort internalVertices(T);
  L := sort leaves(graph(graph(T)));
  n := #L;
  rho := (delete(null, apply(int, i -> if degreeIn(T, i) == 0 then i)))_0;

  -- make source and target rings
  KK := opts.CoefficientRing;
  R := ssmFourierParameterRing(T, CoefficientRing => KK);
  S := if opts.SourceRing === null then KK[apply(ssmStates(n), i -> q_i)] else opts.SourceRing;

  phi := for leafState in ssmStates(n) list(

    upperInd := apply(leafState, k -> k_0);
    lowerInd := apply(leafState, k -> k_1);
    
    sum for intState in (toList(#int:0)..toList(#int:1)) list(

      state := hashTable(apply(L, i -> i => lowerInd_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));

      (r_(state#rho))*(edges(T) / (e ->  (a_(e, sum(apply(leafDescendants(e_1, L, T), l -> upperInd_(l-1))) % 2, state#(e_0), state#(e_1)))) // product)
      )
    );

  if not opts.UseStochasticParameters then return map(R, S, phi);

  subRules := {r_0 => 1}|(flatten for e in edges(T) list {a_(e, 0, 0, 1) => 1 - a_(e, 0, 0, 0), a_(e, 0, 1, 1) => 1 - a_(e, 0, 1, 0)});
  phi = apply(phi, i -> sub(i, subRules));

  return map(R, S, phi)
  );


ssmFourierConditional = (reticulationVertex, N, R) -> (

  w := first toList children(N, reticulationVertex);
  (u, v) := toSequence sort toList parents(N, reticulationVertex);
  states := {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

  condDist = new MutableHashTable;

  for x in states do for j in states do for k in states do(

    if (x_0 + j_0 + k_0) % 2 != 0 then(
      condDist#(x, j, k) = 0;
      continue;
    );
    
    if x_0 == 0 and j_0 == 1 and k_0 == 1 then(
      condDist#(x, j, k) = 0;
      continue;
    );

    if x_0 == 0 and (j_0 + k_0) % 2 == 0 then condDist#(x, j, k) = a_({u, w}, 0, j_1, x_1) + a_({v, w}, 0, k_1, x_1);

    if x_0 == 1 and (j_0 + k_0) % 2 == 1 then(

      if j_0 == 1 then condDist#(x, j, k) = a_({u, w}, 1, j_1, x_1);
      if k_0 == 1 then condDist#(x, j, k) = a_({v, w}, 1, k_1, x_1);
    );
  );
  
  return condDist
)


n = 3
kappa = 4
T = digraph({1,2,3, 4}, {{3, 4}, {2, 4}, {4, 1}})
R = ssmFourierParameterRing({{3, 4}, {2, 4}}, T)
ssmFourierConditional(4, T, R)

-- N, a digraph representing a phylogenetic network
-- creates the parametrization of the strand symmetric model on N in the fourier coordinates
-- By default the stochastic restrictions are enforced but the parametrization is re-homogenized 
ssmFourierNetworkParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true, CoefficientRing => QQ})
ssmFourierNetworkParametrization (List, Digraph) := RingMap => opts -> (reticulationEdges, N) -> (

  -- compute internal vertices, leaves, and displayed trees
	int := sort internalVertices(N);
	L := sort delete(null, vertices(N) / (i -> if degree(N, i) == 1 then i));
	n := #L;
	retVert := (reticulationEdges / set // intersect // toList )_0;
	retLeaf := (toList children(N, retVert))_0;
	(u, v) := toSequence sort toList parents(N, retVert);
  T := deleteEdges(N, reticulationEdges);
  rho := (delete(null, apply(int, i -> if degreeIn(N, i) == 0 then i)))_0;
  states := {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

  -- make source and target rings
  KK := opts.CoefficientRing;
  R = ssmFourierParameterRing(reticulationEdges, N, CoefficientRing => KK);
  S := if opts.SourceRing === null then KK[apply(ssmStates(n), i -> q_i)] else opts.SourceRing;

  -- make the conditional distribution
  condDist := ssmFourierConditional(retVert, N, R);


	phi := for leafState in ssmStates(n) list(

      -- compute the induced state at each internal vertex
	    upperInd := new MutableHashTable from apply(L, leafState, (l, k) -> {l, k_0});
      for i in int do upperInd#i = sum(apply(leafDescendants(i, L, T), l -> upperInd#l)) % 2;
	    lowerInd := apply(leafState, k -> k_1);
	    
	    sum for intState in (toList(#int:0)..toList(#int:1)) list(

	      state := hashTable(apply(L, i -> i => lowerInd_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));

	      phiT := edges(T) / (e -> (a_(e, upperInd#(e_1), state#(e_0), state#(e_1)))) // product;
	      (r_(state#rho))*phiT*(sum flatten for j in states list for k in states list condDist#(leafState_retLeaf, {j_0, state#u}, {k_0, state#v}))
	      )
    	);

	if not opts.UseStochasticParameters then return map(R, S, phi);

  subRules := {r_0 => 1}|(flatten for e in edges(N) list if e == {5, 1} then continue else {a_(e, 0, 0, 1) => 1 - a_(e, 0, 0, 0), a_(e, 0, 1, 1) => 1 - a_(e, 0, 1, 0)});
  subRules = subRules|contractRetLeaf;
  phi = apply(phi, i -> sub(i, subRules));

  return map(R, S, phi)
	);


-- A,B, a pair of list which form a bipartition of [n] representing a split
-- S, the ring of the Strand-Symmetric model on a n-leaf network
-- returns the flattening matrix corresponding to the partition A|B
flat = (A, B, S) -> (

    n := #(last baseName S_0);
    indA := ssmStates(#A);
    indB := ssmStates(#B);
    x := first baseName S_0;
    
    return matrix for a in indA list(

        for b in indB list(

            curInd := new MutableHashTable from apply(#A, i -> A_i => a_i) | apply(#B, i -> B_i => b_i);

            (x_(apply(n, i -> curInd#(i+1))))_S
            )
        );
    );

