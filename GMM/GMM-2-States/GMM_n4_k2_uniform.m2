load("../../PhylogeneticNetworks.m2")
load("../../MatroidDistinguish.m2")
load("../GMM_Networks.m2")


-- create the ring and the parametrization of the 2-state GMM on a 4-sunlet with reticulation 1 and leaves labelled cyclically
n = 4
k = 2
S = pRing(k, n)
N = digraph(toList(1..8), {{5,1}, {6,2}, {7,3}, {8,4}, {6,5}, {8,5}, {7,8}, {7,6}});
phi = gmmNetworkReParametrization(k, {{6, 5}, {8, 5}}, N, SourceRing => S);

J = jacobian matrix phi;
nJ = specialize J

-- check that every set of size 14 is a basis
all(subsets(16, 14), A -> rank nJ_A == 14)
