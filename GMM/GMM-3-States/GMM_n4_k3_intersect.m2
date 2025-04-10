load("../../PhylogeneticNetworks.m2")
load("../../MatroidDistinguish.m2")
load("../GMM_Networks.m2")
needsPackage "NumericalImplicitization")


-- create the ring for 2-state GMM on a 4-leaf network
n = 4
k = 3
S = pRing(k, n)


-- create the network and its parametrization and jacobian
N = digraph(toList(1..8), {{5,1}, {6,2}, {7,3}, {8,4}, {6,5}, {8,5}, {7,8}, {7,6}});
phi = gmmNetworkReParametrization(k, {{6, 5}, {8, 5}}, N, SourceRing => S);
jac = specialize jacobian matrix phi;

-- create the two ideals of minors corresponding to the two displayed trees of a 4-leaf sunlet
flat12 = flat({1, 2}, {3, 4}, S)
flat14 = flat({1, 4}, {2, 3}, S)
cutoff = 5;
flat12' = submatrix(flat12, 0..cutoff, 0..cutoff)
flat14' = submatrix(flat14, 0..cutoff, 0..cutoff)
I1 = minors(k+1, flat12');
I2 = minors(k+1, flat14');

-- intersect these ideals
gbTrace = 3
intersect(I1, I2)


supp = unique(support(submatrix(flat12, toList(0..5), toList(0..5))) | support(submatrix(flat14, toList(0..5), toList(0..5))));
rank (jac_(supp / index))