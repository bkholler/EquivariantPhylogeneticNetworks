load("../../PhylogeneticNetworks.m2")
load("../GMM_Networks.m2")


-- create the ring corresponding to the 2-state general Markov model on a 4-leaf network
k = 2
n = 4
S = pRing(k, n)


-- make all 4 cycles and their corresponding parameterizations
fourCycles = generate4Cycles()
fourCycleMaps = hashTable apply(fourCycles, N -> {N, gmmNetworkParametrization(k, {{6, 5}, {8, 5}}, N,  SourceRing => S)});


-- make our polynomial f = alpha(g)
-- we will show this polynomial belongs to the kernel of the parameterization of a network N if and only if the nontrivial splits of the displayed trees of N are 12|34 and 14|23
flat12 = flat({1,2}, {3,4}, S);
flat14 = flat({1,4}, {2,3}, S);
inds = {{({0, 1, 2},{0, 1, 2}), ({1, 2, 3},{1, 2, 3})}, {({0, 1, 2},{0, 1, 3}), ({0, 2, 3},{1, 2, 3})}, {({0, 1, 2},{0, 2, 3}), ({1, 2, 3},{0, 2, 3})}, {({0, 1, 2},{1, 2, 3}), ({0, 2, 3},{0, 2, 3})}, {({0, 2, 3},{0, 1, 2}), ({0, 1, 3},{1, 2, 3})}, {({0, 2, 3},{0, 1, 3}), ({0, 1, 2},{1, 2, 3})}, {({0, 2, 3},{0, 2, 3}), ({0, 1, 3},{0, 2, 3})}, {({0, 2, 3},{1, 2, 3}), ({0, 1, 2},{0, 2, 3})}};
coeffs = {1, -1, 1, -1, -1, 1, -1, 1};
f = fillTemplate(inds, coeffs, flat12, flat14);


-- first, we quickly find networks such that f does not vanish on the model
-- we do this by sampling points from the model using the parameterization and simply checking that they do not evaluate to 0
-- every network in the following list cannot have f in its vanishing ideal but, since this is probabilistic the converse is only true with high probability
nonvanishingNetworks = time select(fourCycles, N -> not probabilisticMembershipTest(10000, f, fourCycleMaps#N))


