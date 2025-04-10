load("../../PhylogeneticNetworks.m2")
load("../../MatroidDistinguish.m2")
load("../GMM_Networks.m2")
needsPackage "MultigradedImplicitization"


-- create the ring for 2-state GMM on a 4-leaf network
n = 4
k = 3
S = pRing(k, n)


-- create the network and its parametrization and jacobian
fourSunlets = generate4Cycles()
N1 = fourSunlets_0
phi1 = gmmNetworkReParametrization(k, {{6, 5}, {8, 5}}, N1, SourceRing => S);
jac1 = specialize jacobian matrix phi1;

N2 = fourSunlets_1
phi2 = gmmNetworkReParametrization(k, {{6, 5}, {8, 5}}, N2, SourceRing => S);
jac2 = specialize jacobian matrix phi2;

-- create the two ideals of minors corresponding to the two displayed trees of the sunlet N1
flat12 = flat({1, 2}, {3, 4}, S)
flat14 = flat({1, 4}, {2, 3}, S)
cutoff = 5;
flat12' = submatrix(flat12, 0..cutoff, 0..cutoff)
flat14' = submatrix(flat14, 0..cutoff, 0..cutoff)
I1 = minors(k+1, flat12');
I2 = minors(k+1, flat14');

-- the support of these two small flattenings form a dependent set of our matroid since it has rank 42
-- however it is also dependent in the matroid for the network N2
supp1214 = unique(support(flat12')|support(flat14'))
rank (jac1_(supp1214 / index))
rank (jac2_(supp1214 / index))

-- this dependence is most likely due to the fact that the we can get the same support if we take flat13 instead of flat14
flat13 = flat({1, 3}, {2, 4}, S)
flat13' = submatrix(flat13, 0..cutoff, 0..cutoff)
supp1213 = unique(support(flat12')|support(flat13'))
rank (jac1_(supp1213 / index))
rank (jac2_(supp1213 / index))
sort(supp1213) == sort(supp1214)


-- we can try to compute the polynomials in these 48 variables and see if
C = QQ[supp1214]
psi1 = map(target phi1, C, (matrix phi1)_(supp1214 / index));
A = transpose matrix for x in gens(S) list(

	l = last baseName x; 

	flatten for l' in l list apply(k, i -> if i == l' then 1 else 0)
);
B' = submatrix'(A, toList(0..k-1),)
B = B'_(supp1214/index)

componentsOfKernel(3, psi1, Grading => B, UseInterpolation => true, ReduceFirst => true)