load("../../PhylogeneticNetworks.m2")
load("../GMM_Networks.m2")


-- make the network
n = 3
k = 4
N = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 2}}, {4, {1}}})


-- make the parametrization of 4-state GMM and its Jacobian
phi = gmmNetworkReParametrization(k, {{6, 4}, {5, 4}}, N);
J = jacobian matrix phi; 

-- compute the dimension over a finite field
-- since it achieves the expected dimension at a random point, it must have the expected dimension
KK = ZZ/nextPrime(10000000);
evalJ = sub(J, apply(gens ring J, i -> i => random(KK)));

-- now we show that this is the uniform matroid
-- to do this we simply check that every subset of the columns of size 61 has full rank
dependentSets = for B in subsets(64, 61) list if rank(evalJ_B) < 61 then B;
