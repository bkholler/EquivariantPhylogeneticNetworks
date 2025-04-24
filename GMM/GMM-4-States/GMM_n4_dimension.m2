load("../../PhylogeneticNetworks.m2")
load("../GMM_Networks.m2")


-- make the network
n = 4
k = 4
N = (generate4Cycles())_0


-- make the parametrization of 4-state GMM and its Jacobian
-- the following code takes about 5 minutes
phi = time gmmNetworkReParametrization(k, {{8, 5}, {6, 5}}, N);
J = jacobian matrix phi; 


-- compute the dimension over a finite field
-- since it achieves the expected dimension at a random point, it must have the expected dimension
KK = ZZ/nextPrime(10000000);
evalJ = sub(J, apply(gens ring J, i -> i => random(KK)));
rank evalJ

