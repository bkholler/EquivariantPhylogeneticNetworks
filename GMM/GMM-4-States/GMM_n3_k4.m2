load("../../PhylogeneticNetworks.m2")
load("../../MatroidDistinguish.m2")
load("../GMM_Networks.m2")


-- create the ring and the parametrization of the 4-state GMM on a 3-sunlet with reticulation 1 and leaves labelled cyclically
n = 3
k = 4
N = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 2}}, {4, {1}}})
S = pRing(k, n)

-- make the parametrization of 4-state GMM and its Jacobian
phi = gmmNetworkReParametrization(k, {{6, 4}, {5, 4}}, N);
J = jacobian matrix phi; 
evalJ = specialize(J, CoefficientRing => ZZ/nextPrime(2^31));

-- compute the dimension over a finite field
-- since it achieves the expected dimension at a random point, it must have the expected dimension
rank evalJ

-- check that every set of size 61 is a basis
-- this shows that the algebraic matroid of the 3-sunlet network is the uniform matroid
-- thus the matroid-based approach cannot be used to distinguish these networks
all(subsets(64, 61), A -> rank evalJ_A == 61)