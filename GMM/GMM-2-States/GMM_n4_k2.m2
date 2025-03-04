load("../../PhylogeneticNetworks.m2")
load("../GMM_Networks.m2")
needsPackage "Permutations"


-- create the ring corresponding to the 2-state general Markov model on a 4-leaf network
k = 2
n = 4
S = pRing(k, n)


-- make all 4-sunlets by permuting the leaves
fourSunlets = generate4Cycles()


-- load the ideal of the cyclically labelled 4-sunlet with reticulation at 1. This is the network fourSunlets_0
load("GMM_n4_k2_ideal.m2")


-- make a hash table where the keys are the four sunlets and the values are their ideal, created by permuting the indices of the original ideal
fourSunletIdeals = hashTable for i in generate4CyclesWithPerms() list(

    (sigma, N) := i;
    tau := (sigma / (i -> i + 1)) // permutation // inverse // toList / (i -> i -1);
    permuteVars :=  apply(gens S, j -> j => p_(toSequence (last baseName j)_tau));
    J := sub(I, permuteVars);
    {N, J}
);


-- To verify that all 4-sunlets are identifiable from each other, we just need to show that every pair of ideals is distinct
all(subsets(keys fourSunletIdeals, 2), netPair ->  fourSunletIdeals#(netPair_0) != fourSunletIdeals#(netPair_1))

