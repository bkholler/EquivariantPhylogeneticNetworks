load("../../PhylogeneticNetworks.m2")
load("../GMM_Networks.m2")

-- create the ring and the vanishing ideal of the 2-state GMM on a 4-sunlet with reticulation 1 and leaves labelled cyclically
n = 4
k = 2
S = pRing(k, n)
I = ideal(p_(0,0,1,1)*p_(0,1,0,0)*p_(0,1,0,1)*p_(1,0,0,0)-p_(0,0,1,0)*p_(0,1,0,1)^2*p_(1,0,0,0)+p_(0,0,0,1)*p_(0,1,0,1)*p_(0,1,1,0)*p_(1,0,0,0)-p_(0,0,0,1)*p_(0,1,0,0)*p_(0,1,1,1)*p_(1,0,0,0)-p_(0,0,1,1)*p_(0,1,0,0)^2*p_(1,0,0,1)+p_(0,0,1,0)*p_(0,1,0,0)*p_(0,1,0,1)*p_(1,0,0,1)-p_(0,0,0,0)*p_(0,1,0,1)*p_(0,1,1,0)*p_(1,0,0,1)+p_(0,0,0,0)*p_(0,1,0,0)*p_(0,1,1,1)*p_(1,0,0,1)-p_(0,0,0,1)*p_(0,1,0,0)*p_(0,1,0,1)*p_(1,0,1,0)+p_(0,0,0,0)*p_(0,1,0,1)^2*p_(1,0,1,0)+p_(0,0,0,1)*p_(0,1,0,0)^2*p_(1,0,1,1)-p_(0,0,0,0)*p_(0,1,0,0)*p_(0,1,0,1)*p_(1,0,1,1)+p_(0,0,0,1)*p_(0,0,1,0)*p_(0,1,0,1)*p_(1,1,0,0)-p_(0,0,0,0)*p_(0,0,1,1)*p_(0,1,0,1)*p_(1,1,0,0)-p_(0,0,0,1)^2*p_(0,1,1,0)*p_(1,1,0,0)+p_(0,0,0,0)*p_(0,0,0,1)*p_(0,1,1,1)*p_(1,1,0,0)-2*p_(0,0,1,1)*p_(0,1,0,0)*p_(1,0,0,1)*p_(1,1,0,0)+2*p_(0,0,1,0)*p_(0,1,0,1)*p_(1,0,0,1)*p_(1,1,0,0)-2*p_(0,0,0,1)*p_(0,1,1,0)*p_(1,0,0,1)*p_(1,1,0,0)+2*p_(0,0,0,0)*p_(0,1,1,1)*p_(1,0,0,1)*p_(1,1,0,0)+p_(0,1,1,1)*p_(1,0,0,0)*p_(1,0,0,1)*p_(1,1,0,0)-p_(0,1,1,0)*p_(1,0,0,1)^2*p_(1,1,0,0)+p_(0,1,0,1)*p_(1,0,0,1)*p_(1,0,1,0)*p_(1,1,0,0)+2*p_(0,0,0,1)*p_(0,1,0,0)*p_(1,0,1,1)*p_(1,1,0,0)-2*p_(0,0,0,0)*p_(0,1,0,1)*p_(1,0,1,1)*p_(1,1,0,0)-p_(0,1,0,1)*p_(1,0,0,0)*p_(1,0,1,1)*p_(1,1,0,0)-p_(0,0,1,1)*p_(1,0,0,1)*p_(1,1,0,0)^2+p_(0,0,0,1)*p_(1,0,1,1)*p_(1,1,0,0)^2-p_(0,0,0,1)*p_(0,0,1,0)*p_(0,1,0,0)*p_(1,1,0,1)+p_(0,0,0,0)*p_(0,0,1,1)*p_(0,1,0,0)*p_(1,1,0,1)+p_(0,0,0,0)*p_(0,0,0,1)*p_(0,1,1,0)*p_(1,1,0,1)-p_(0,0,0,0)^2*p_(0,1,1,1)*p_(1,1,0,1)+2*p_(0,0,1,1)*p_(0,1,0,0)*p_(1,0,0,0)*p_(1,1,0,1)-2*p_(0,0,1,0)*p_(0,1,0,1)*p_(1,0,0,0)*p_(1,1,0,1)+2*p_(0,0,0,1)*p_(0,1,1,0)*p_(1,0,0,0)*p_(1,1,0,1)-2*p_(0,0,0,0)*p_(0,1,1,1)*p_(1,0,0,0)*p_(1,1,0,1)-p_(0,1,1,1)*p_(1,0,0,0)^2*p_(1,1,0,1)+p_(0,1,1,0)*p_(1,0,0,0)*p_(1,0,0,1)*p_(1,1,0,1)-2*p_(0,0,0,1)*p_(0,1,0,0)*p_(1,0,1,0)*p_(1,1,0,1)+2*p_(0,0,0,0)*p_(0,1,0,1)*p_(1,0,1,0)*p_(1,1,0,1)-p_(0,1,0,0)*p_(1,0,0,1)*p_(1,0,1,0)*p_(1,1,0,1)+p_(0,1,0,0)*p_(1,0,0,0)*p_(1,0,1,1)*p_(1,1,0,1)+p_(0,0,1,1)*p_(1,0,0,0)*p_(1,1,0,0)*p_(1,1,0,1)+p_(0,0,1,0)*p_(1,0,0,1)*p_(1,1,0,0)*p_(1,1,0,1)-p_(0,0,0,1)*p_(1,0,1,0)*p_(1,1,0,0)*p_(1,1,0,1)-p_(0,0,0,0)*p_(1,0,1,1)*p_(1,1,0,0)*p_(1,1,0,1)-p_(0,0,1,0)*p_(1,0,0,0)*p_(1,1,0,1)^2+p_(0,0,0,0)*p_(1,0,1,0)*p_(1,1,0,1)^2+p_(0,0,0,1)^2*p_(0,1,0,0)*p_(1,1,1,0)-p_(0,0,0,0)*p_(0,0,0,1)*p_(0,1,0,1)*p_(1,1,1,0)+2*p_(0,0,0,1)*p_(0,1,0,0)*p_(1,0,0,1)*p_(1,1,1,0)-2*p_(0,0,0,0)*p_(0,1,0,1)*p_(1,0,0,1)*p_(1,1,1,0)-p_(0,1,0,1)*p_(1,0,0,0)*p_(1,0,0,1)*p_(1,1,1,0)+p_(0,1,0,0)*p_(1,0,0,1)^2*p_(1,1,1,0)+p_(0,0,0,1)*p_(1,0,0,0)*p_(1,1,0,1)*p_(1,1,1,0)-p_(0,0,0,0)*p_(1,0,0,1)*p_(1,1,0,1)*p_(1,1,1,0)-p_(0,0,0,0)*p_(0,0,0,1)*p_(0,1,0,0)*p_(1,1,1,1)+p_(0,0,0,0)^2*p_(0,1,0,1)*p_(1,1,1,1)-2*p_(0,0,0,1)*p_(0,1,0,0)*p_(1,0,0,0)*p_(1,1,1,1)+2*p_(0,0,0,0)*p_(0,1,0,1)*p_(1,0,0,0)*p_(1,1,1,1)+p_(0,1,0,1)*p_(1,0,0,0)^2*p_(1,1,1,1)-p_(0,1,0,0)*p_(1,0,0,0)*p_(1,0,0,1)*p_(1,1,1,1)-p_(0,0,0,1)*p_(1,0,0,0)*p_(1,1,0,0)*p_(1,1,1,1)+p_(0,0,0,0)*p_(1,0,0,1)*p_(1,1,0,0)*p_(1,1,1,1),p_(0,0,1,1)*p_(0,1,0,1)*p_(0,1,1,0)*p_(1,0,0,0)-p_(0,0,1,0)*p_(0,1,0,1)*p_(0,1,1,1)*p_(1,0,0,0)-p_(0,0,1,1)*p_(0,1,0,0)*p_(0,1,1,0)*p_(1,0,0,1)+p_(0,0,1,0)*p_(0,1,0,0)*p_(0,1,1,1)*p_(1,0,0,1)-p_(0,0,0,1)*p_(0,1,0,0)*p_(0,1,1,1)*p_(1,0,1,0)+p_(0,0,0,0)*p_(0,1,0,1)*p_(0,1,1,1)*p_(1,0,1,0)+p_(0,0,0,1)*p_(0,1,0,0)*p_(0,1,1,0)*p_(1,0,1,1)-p_(0,0,0,0)*p_(0,1,0,1)*p_(0,1,1,0)*p_(1,0,1,1)-p_(0,0,0,1)*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,1,0,0)+p_(0,0,0,1)*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,1,0,0)-2*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,0,0,1)*p_(1,1,0,0)+2*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,0,0,1)*p_(1,1,0,0)+p_(0,1,1,1)*p_(1,0,0,1)*p_(1,0,1,0)*p_(1,1,0,0)-p_(0,1,1,0)*p_(1,0,0,1)*p_(1,0,1,1)*p_(1,1,0,0)+p_(0,0,0,0)*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,1,0,1)-p_(0,0,0,0)*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,1,0,1)+2*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,0,0,0)*p_(1,1,0,1)-2*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,0,0,0)*p_(1,1,0,1)-p_(0,1,1,1)*p_(1,0,0,0)*p_(1,0,1,0)*p_(1,1,0,1)+p_(0,1,1,0)*p_(1,0,0,0)*p_(1,0,1,1)*p_(1,1,0,1)+p_(0,0,0,1)*p_(0,0,1,1)*p_(0,1,0,0)*p_(1,1,1,0)-p_(0,0,0,0)*p_(0,0,1,1)*p_(0,1,0,1)*p_(1,1,1,0)+2*p_(0,0,0,1)*p_(0,1,0,0)*p_(1,0,1,1)*p_(1,1,1,0)-2*p_(0,0,0,0)*p_(0,1,0,1)*p_(1,0,1,1)*p_(1,1,1,0)-p_(0,1,0,1)*p_(1,0,0,0)*p_(1,0,1,1)*p_(1,1,1,0)+p_(0,1,0,0)*p_(1,0,0,1)*p_(1,0,1,1)*p_(1,1,1,0)-p_(0,0,1,1)*p_(1,0,0,1)*p_(1,1,0,0)*p_(1,1,1,0)+p_(0,0,0,1)*p_(1,0,1,1)*p_(1,1,0,0)*p_(1,1,1,0)+p_(0,0,1,1)*p_(1,0,0,0)*p_(1,1,0,1)*p_(1,1,1,0)-p_(0,0,0,0)*p_(1,0,1,1)*p_(1,1,0,1)*p_(1,1,1,0)-p_(0,0,0,1)*p_(0,0,1,0)*p_(0,1,0,0)*p_(1,1,1,1)+p_(0,0,0,0)*p_(0,0,1,0)*p_(0,1,0,1)*p_(1,1,1,1)-2*p_(0,0,0,1)*p_(0,1,0,0)*p_(1,0,1,0)*p_(1,1,1,1)+2*p_(0,0,0,0)*p_(0,1,0,1)*p_(1,0,1,0)*p_(1,1,1,1)+p_(0,1,0,1)*p_(1,0,0,0)*p_(1,0,1,0)*p_(1,1,1,1)-p_(0,1,0,0)*p_(1,0,0,1)*p_(1,0,1,0)*p_(1,1,1,1)+p_(0,0,1,0)*p_(1,0,0,1)*p_(1,1,0,0)*p_(1,1,1,1)-p_(0,0,0,1)*p_(1,0,1,0)*p_(1,1,0,0)*p_(1,1,1,1)-p_(0,0,1,0)*p_(1,0,0,0)*p_(1,1,0,1)*p_(1,1,1,1)+p_(0,0,0,0)*p_(1,0,1,0)*p_(1,1,0,1)*p_(1,1,1,1),p_(0,0,1,1)*p_(0,1,1,0)*p_(0,1,1,1)*p_(1,0,0,0)-p_(0,0,1,0)*p_(0,1,1,1)^2*p_(1,0,0,0)-p_(0,0,1,1)*p_(0,1,1,0)^2*p_(1,0,0,1)+p_(0,0,1,0)*p_(0,1,1,0)*p_(0,1,1,1)*p_(1,0,0,1)+p_(0,0,1,1)*p_(0,1,0,1)*p_(0,1,1,0)*p_(1,0,1,0)-p_(0,0,1,1)*p_(0,1,0,0)*p_(0,1,1,1)*p_(1,0,1,0)-p_(0,0,0,1)*p_(0,1,1,0)*p_(0,1,1,1)*p_(1,0,1,0)+p_(0,0,0,0)*p_(0,1,1,1)^2*p_(1,0,1,0)-p_(0,0,1,0)*p_(0,1,0,1)*p_(0,1,1,0)*p_(1,0,1,1)+p_(0,0,0,1)*p_(0,1,1,0)^2*p_(1,0,1,1)+p_(0,0,1,0)*p_(0,1,0,0)*p_(0,1,1,1)*p_(1,0,1,1)-p_(0,0,0,0)*p_(0,1,1,0)*p_(0,1,1,1)*p_(1,0,1,1)-p_(0,0,1,1)^2*p_(0,1,1,0)*p_(1,1,0,0)+p_(0,0,1,0)*p_(0,0,1,1)*p_(0,1,1,1)*p_(1,1,0,0)-2*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,0,1,1)*p_(1,1,0,0)+2*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,0,1,1)*p_(1,1,0,0)+p_(0,1,1,1)*p_(1,0,1,0)*p_(1,0,1,1)*p_(1,1,0,0)-p_(0,1,1,0)*p_(1,0,1,1)^2*p_(1,1,0,0)+p_(0,0,1,0)*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,1,0,1)-p_(0,0,1,0)^2*p_(0,1,1,1)*p_(1,1,0,1)+2*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,0,1,0)*p_(1,1,0,1)-2*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,0,1,0)*p_(1,1,0,1)-p_(0,1,1,1)*p_(1,0,1,0)^2*p_(1,1,0,1)+p_(0,1,1,0)*p_(1,0,1,0)*p_(1,0,1,1)*p_(1,1,0,1)+p_(0,0,1,1)^2*p_(0,1,0,0)*p_(1,1,1,0)-p_(0,0,1,0)*p_(0,0,1,1)*p_(0,1,0,1)*p_(1,1,1,0)+p_(0,0,0,1)*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,1,1,0)-p_(0,0,0,0)*p_(0,0,1,1)*p_(0,1,1,1)*p_(1,1,1,0)-2*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,0,0,1)*p_(1,1,1,0)+2*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,0,0,1)*p_(1,1,1,0)+p_(0,1,1,1)*p_(1,0,0,1)*p_(1,0,1,0)*p_(1,1,1,0)+2*p_(0,0,1,1)*p_(0,1,0,0)*p_(1,0,1,1)*p_(1,1,1,0)-2*p_(0,0,1,0)*p_(0,1,0,1)*p_(1,0,1,1)*p_(1,1,1,0)+2*p_(0,0,0,1)*p_(0,1,1,0)*p_(1,0,1,1)*p_(1,1,1,0)-2*p_(0,0,0,0)*p_(0,1,1,1)*p_(1,0,1,1)*p_(1,1,1,0)-p_(0,1,1,1)*p_(1,0,0,0)*p_(1,0,1,1)*p_(1,1,1,0)-p_(0,1,0,1)*p_(1,0,1,0)*p_(1,0,1,1)*p_(1,1,1,0)+p_(0,1,0,0)*p_(1,0,1,1)^2*p_(1,1,1,0)+p_(0,0,1,1)*p_(1,0,1,0)*p_(1,1,0,1)*p_(1,1,1,0)-p_(0,0,1,0)*p_(1,0,1,1)*p_(1,1,0,1)*p_(1,1,1,0)-p_(0,0,1,1)*p_(1,0,0,1)*p_(1,1,1,0)^2+p_(0,0,0,1)*p_(1,0,1,1)*p_(1,1,1,0)^2-p_(0,0,1,0)*p_(0,0,1,1)*p_(0,1,0,0)*p_(1,1,1,1)+p_(0,0,1,0)^2*p_(0,1,0,1)*p_(1,1,1,1)-p_(0,0,0,1)*p_(0,0,1,0)*p_(0,1,1,0)*p_(1,1,1,1)+p_(0,0,0,0)*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,1,1,1)+2*p_(0,0,1,1)*p_(0,1,1,0)*p_(1,0,0,0)*p_(1,1,1,1)-2*p_(0,0,1,0)*p_(0,1,1,1)*p_(1,0,0,0)*p_(1,1,1,1)-2*p_(0,0,1,1)*p_(0,1,0,0)*p_(1,0,1,0)*p_(1,1,1,1)+2*p_(0,0,1,0)*p_(0,1,0,1)*p_(1,0,1,0)*p_(1,1,1,1)-2*p_(0,0,0,1)*p_(0,1,1,0)*p_(1,0,1,0)*p_(1,1,1,1)+2*p_(0,0,0,0)*p_(0,1,1,1)*p_(1,0,1,0)*p_(1,1,1,1)-p_(0,1,1,0)*p_(1,0,0,1)*p_(1,0,1,0)*p_(1,1,1,1)+p_(0,1,0,1)*p_(1,0,1,0)^2*p_(1,1,1,1)+p_(0,1,1,0)*p_(1,0,0,0)*p_(1,0,1,1)*p_(1,1,1,1)-p_(0,1,0,0)*p_(1,0,1,0)*p_(1,0,1,1)*p_(1,1,1,1)-p_(0,0,1,1)*p_(1,0,1,0)*p_(1,1,0,0)*p_(1,1,1,1)+p_(0,0,1,0)*p_(1,0,1,1)*p_(1,1,0,0)*p_(1,1,1,1)+p_(0,0,1,1)*p_(1,0,0,0)*p_(1,1,1,0)*p_(1,1,1,1)+p_(0,0,1,0)*p_(1,0,0,1)*p_(1,1,1,0)*p_(1,1,1,1)-p_(0,0,0,1)*p_(1,0,1,0)*p_(1,1,1,0)*p_(1,1,1,1)-p_(0,0,0,0)*p_(1,0,1,1)*p_(1,1,1,0)*p_(1,1,1,1)-p_(0,0,1,0)*p_(1,0,0,0)*p_(1,1,1,1)^2+p_(0,0,0,0)*p_(1,0,1,0)*p_(1,1,1,1)^2);


-- verify that I is prime and has the right dimension
isPrime I
dim(I) == 14


-- setup ring of variables corresponding to minors of flattenings
R = QQ[(flatten for A in subsets(n, 2) list for B in subsets(n, 2) list x_(A,B)) | (flatten for A in subsets(n, 2) list for B in subsets(n, 2) list y_(A,B))]
flat12 = flat({1, 2}, {3, 4}, S)
flat14 = flat({1, 4}, {2, 3}, S)
images = (flatten for A in subsets(n, 2) list for B in subsets(n, 2) list det(flat12_B^A))|(flatten for A in subsets(n, 2) list for B in subsets(n, 2) list det(flat14_B^A));
psi = map(S, R, images);

-- Create the matrix of minors shown in Figure 2
A = matrix {{-y_({0, 1},{1, 3})+y_({1, 2},{1, 3})-y_({0, 3},{1, 3})-y_({2,3},{1, 3}), y_({0, 1},{1, 3})-y_({2, 3},{1, 3})}, 
{-y_({0, 1},{1, 2})-y_({0,1},{0, 3})+y_({1, 2},{1, 2})+y_({1, 2},{0, 3})-y_({0, 3},{1, 2})-y_({0, 3},{0,3})-y_({2, 3},{1, 2})-y_({2, 3},{0, 3}), y_({0, 1},{1, 2})+y_({0, 1},{0,3})-y_({2, 3},{1, 2})-y_({2, 3},{0, 3})}, 
{-y_({0, 1},{0, 2})+y_({1, 2},{0,2})-y_({0, 3},{0, 2})-y_({2, 3},{0, 2}), y_({0, 1},{0, 2})-y_({2, 3},{0, 2})}}

-- Verify that the 2x2 minors of this matrix are the generators of I after applying psi to the entries of A
I == minors(2, psi(A))


-- verify that the polynomials in I vanish on the 2-state GMM model
-- First, we apply phi and psi to A to get a matrix with entries in the parameter space of the model. We then verify that it has rank 1, thus all of the 2x2 minors are indeed in the kernel
N = digraph(toList(1..8), {{5,1}, {6,2}, {7,3}, {8,4}, {6,5}, {8,5}, {7,8}, {7,6}});
phi = gmmNetworkReParametrization(k, {{6, 5}, {8, 5}}, N, SourceRing => S);

-- The code below takes about 20 seconds
rank(phi(psi(A)))