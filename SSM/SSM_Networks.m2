load "../PhylogeneticNetworks.m2"


-- generate all leaf-labellings for the SSM model on n leaves
ssmStates = n -> (

  states := toList(toList(n:{0,0})..toList(n:{1,1}));

  return delete(null, apply(states, state -> if sum(apply(state, ind -> ind_0)) % 2 == 0 then state))
  )

-- N, a digraph representing a phylogenetic network
-- creates the ring of fourier parameters for the SSM model on N
ssmParameterRing = method(Options => {aVariableName => "a"})
ssmParameterRing Digraph := Ring => opts -> N -> (

  L := sort delete(null, apply(vertices(N), i -> if degree(N, i) == 1 then i));
  n := #L;
  a := getSymbol opts.aVariableName;

  return QQ[s, t, flatten for e in edges(T) list a_(e, 0, 0, 0)..a_(e, 1, 1, 1)];
  )


-- T, a digraph representing a phylogenetic tree
-- creates the parametrization of the strand symmetric model on T in the fourier coordinates
ssmTreeParametrization = method(Options => {SourceRing => null, UseStochasticParameters => true})
ssmTreeParametrization Digraph := RingMap => opts -> T -> (

  int := sort internalVertices(T);
  L := sort leaves(graph(graph(T)));
  n := #L;
  R := ssmParameterRing(T);
  S := if opts.SourceRing === null then QQ[apply(ssmStates(n), i -> q_i)] else opts.SourceRing;

  phi := for leafState in ssmStates(n) list(

    upperInd := apply(leafState, k -> k_0);
    lowerInd := apply(leafState, k -> k_1);
    
    sum for intState in (toList(#int:0)..toList(#int:1)) list(

      state := hashTable(apply(L, i -> i => lowerInd_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));

      edges(T) / (e ->  (a_(e, sum(apply(leafDescendants(e_1, L, T), l -> upperInd_(l-1))) % 2, state#(e_0), state#(e_1)))) // product
      )
    );

  if not opts.UseStochasticParameters then return map(R, S, phi);

  subRules := flatten for e in edges(T) list {a_(e, 0, 0, 1) => 1 - a_(e, 0, 0, 0)};
  phi = apply(phi, i -> s*sub(i, subRules));

  return map(R, S, phi)
  );

ssmNetParam = (N, retEdges) -> (

	int := sort internalVertices(N);
  L := sort delete(null, apply(vertices(N), i -> if degree(N, i) == 1 then i));
  n := #L;
  S := ssmParamRing(N, "a");
  m := getSymbol "a";
  T1 := deleteEdges(N, {retEdges_0});
  T2 := deleteEdges(N, {retEdges_1});

	phi := for leafState in ssmStates(n) list(

	    upperInd := apply(leafState, k -> k_0);
	    lowerInd := apply(leafState, k -> k_1);
	    
	    sum for intState in (toList(#int:0)..toList(#int:1)) list(

	      state := hashTable(apply(L, i -> i => lowerInd_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));


	      phi1 := edges(T1) / (e -> try (m_(e, sum(apply(leafDescendants(T1, L, e_1), l -> upperInd_(l-1))) % 2, state#(e_1), state#(e_0)))_S) // product;
	      phi2 := edges(T2) / (e -> try (m_(e, sum(apply(leafDescendants(T2, L, e_1), l -> upperInd_(l-1))) % 2, state#(e_1), state#(e_0)))_S) // product;

	      (phi1+phi2)
	      )
  
    	);

	return phi;
	);


flat = (n, A, B, S) -> (

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



interleave = (curInd, pos, inds) -> (

  liftInd := curInd;

  for i from 0 to #pos - 1 do(

    liftInd = insert(pos_i, inds_i, liftInd);
    );

  return liftInd;
  );

liftPoly = (pos, inds, f, oldR, newR) -> (

  oldInds := gens(oldR) / (i -> last baseName i);
  newInds := oldInds / (i -> interleave(i, pos, inds));

  x := first baseName newR_0;
  newVars := apply(newInds, i -> (x_i)_newR);

  subMap := map(newR, oldR, newVars);

  return subMap(f)
  );




