load("../../PhylogeneticNetworks.m2")
load("../../MatroidDistinguish.m2")
load("../GMM_Networks.m2")
needsPackage "MultigradedImplicitization"

N = digraph(toList(1..10), {{9, {10, 8, 4}}, {8, {3, 7}}, {7, {2, 6}}, {10, {5, 6}}, {6, {1}}})
k = 2

-- since the jacobian is full rank, the resulting variety is the whole space
phi = gmmNetworkReParametrization(k, {{10, 6}, {7, 6}}, N, UseStochasticParameters => true);
J = jacobian matrix phi;

-- this computes the rank of the jacobian with random parameter values many times
-- it is at least 19 and thus the model itself has dimension 18. this takes about 30 seconds
tally apply(1000, i -> rank specialize J)


-- we can also compute the ideal of this model over a finite field
-- the resulting polynomials lie in the true ideal with high probability
-- however we cannot verify that they generate the full ideal
KK = ZZ/32003
R = KK[gens target phi]
S = KK[gens source phi]
psi = map(R, S, sub(matrix phi, R))
A = transpose matrix for x in gens(S) list(

	l = last baseName x; 

	flatten for l' in l list apply(k, i -> if i == l' then 1 else 0)
);

B = submatrix'(A, toList(0..k-1),)

-- this takes about 5 minutes
G = componentsOfKernel(4, psi, UseInterpolation => true, Grading => B);
I = ideal delete(null, flatten values G);