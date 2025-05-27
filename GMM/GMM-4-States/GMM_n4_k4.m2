load("../../PhylogeneticNetworks.m2")
load("../../MatroidDistinguish.m2")
load("../GMM_Networks.m2")


-- create the ring corresponding to the 4-state general Markov model on a 4-leaf network
k = 4
n = 4
S = pRing(k, n)


-- make all 4-sunlets by permuting the leaves
-- then select the four sunlets which have reticulation at 1
-- it suffices to distinguish these
fourSunlets = generate4Cycles()
ret1Networks = select(fourSunlets, N -> member({5, 1}, edges(N)))
(N1, N2, N3) = toSequence ret1Networks

-- the following two lines which make the maps and jacobians take about 10-20 minutes to run
-- make a hash table where the keys are the four sunlets and the values are their ideal, created by permuting the indices of the original ideal
fourSunletMaps = hashTable for N in ret1Networks list {N, gmmNetworkReParametrization(k, {{6, 5}, {8, 5}}, N)};

-- To verify that all 4-sunlets are identifiable from each other, we use the matroid-based approach from "Identifiability in Phylogenetics using Algebraic Matroids"
-- However, we do not select random subsets but instead use a particular subset which we have determined already
fourSunletJacs = applyValues(fourSunletMaps, phi -> jacobian matrix phi);

-- compute the dimension over a finite field
-- since it achieves the expected dimension at a random point, it must have the expected dimension
rank specialize fourSunletJacs#N1

-- now we pick a special subset of the variables such that the submatrix of the jacobian of N1 will drop rank 
-- but the corresponding submatrix of the jacobian of both N2 and N3 will be full rank 
flat12 = flat({1, 2}, {3, 4}, S)
C = apply(support submatrix(flat12, {0,1,4,5,8,9, 12, 13}, {0,1,4,5,8,9, 12, 13}), index)

rank specialize (fourSunletJacs#N1)_C
rank specialize (fourSunletJacs#N2)_C
rank specialize (fourSunletJacs#N3)_C

-- we now choose a large finite field according to the degree of the maximal minors of J(phi)
-- we then repeat the rank computation for ret1Networks_0 several times
-- this ensures the probability of identifiability is > 1 - 10^(-30)
alpha = max(flatten(apply(flatten entries fourSunletJacs#N1, degree)));
p = nextPrime((#C)*alpha + 10^10);
KK = ZZ/p;
l = ceiling(log(2, 10^(-30))/log(2, alpha/p));
-- when the following function returns true, then the statement rank (fourSunletJac#N1)_C < #C
-- holds with probability > 1 - 10^(-30)
schwartzZippelCertify(p, l, #C, (fourSunletJacs#N1)_C)


