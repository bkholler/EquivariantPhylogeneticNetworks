load "PhylogeneticNetworks.m2"


-- generate all leaf-labellings for the SSM model on n leaves
ssmStates = n -> (

  states := toList(toList(n:{0,0})..toList(n:{1,1}));

  return delete(null, apply(states, state -> if sum(apply(state, ind -> ind_0)) % 2 == 0 then state))
  )


-- create the ring of parameters for the SSM model on a tree T using symbol paramSym
ssmParamRing = (T, paramSym) -> (

  L := sort delete(null, apply(vertices(T), i -> if degree(T, i) == 1 then i));
  n := #L;
  m := getSymbol paramSym;

  indSet := flatten for e in edges(T) list(


    flatten flatten for j from 0 to 1 list for s from 0 to 1 list for t from 0 to 1 list (toList(e), j, s, t)
    );

  return R := QQ[apply(indSet, i -> m_i)];
  )



ssmTreeParam = (T, paramSym) -> (

  int := sort intVertices(T);
  L := sort leaves(graph(graph(T)));
  n := #L;
  S := ssmParamRing(T, paramSym);
  m := getSymbol paramSym;
  use S;
  
  phi := for leafState in ssmStates(n) list(

    upperInd := apply(leafState, k -> k_0);
    lowerInd := apply(leafState, k -> k_1);
    
    sum for intState in (toList(#int:0)..toList(#int:1)) list(

      state := hashTable(apply(L, i -> i => lowerInd_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));


      edges(T) / (e -> try (m_(e, sum(apply(leafDescendants(T, L, e_1), l -> upperInd_(l-1))) % 2, state#(e_1), state#(e_0)))_S) // product
      )
  
    );

  return phi
  )



ssmNetParam = (N, retEdges) -> (

	int := sort intVertices(N);
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



fillTemplate = (inds, coeffs, flat1, flat2, S) -> (

  return sum for i from 0 to #inds - 1 list(
	
    ind  = inds_i;
    (R1, C1) = ind_0;
    (R2, C2)= ind_1;

    (coeffs_i)*(det(flat1_C1^R1))*(det(flat2_C2^R2))
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




