load("../../PhylogeneticNetworks.m2")
load("../GMM_Networks.m2")


-- create the ring corresponding to the 2-state general Markov model on a 4-leaf network
k = 2
n = 4
S = pRing(k, n)


-- make all 4 cycles and their corresponding parameterizations
fourCycles = generate4Cycles()
fourCycleMaps = hashTable apply(fourCycles, N -> {N, gmmNetworkParametrization(k, {{6, 5}, {8, 5}}, N,  SourceRing => S)});


-- load the ideal of the cyclically labelled 4-cycle with reticulation at 1. This is the network fourCycles_0
load("GMM_n4_k2_ideal.m2")


-- first, we quickly find networks such that f does not vanish on the model
-- we do this by sampling points from the model using the parameterization and simply checking that they do not evaluate to 0
-- if the probabilisticMembershipshipTest(t, f, phi) returns false, then it is unconditionally true that f is not in the kernel of phi. The output is only uncertain if it returns true.
-- As a result, we see that the 11 other distinct 4-cycles are all distinguishable from fourCycles_0
distinguishableNetworks = select(fourCycles, N -> any(I_*, f -> not probabilisticMembershipTest(100, f, fourCycleMaps#N)))
