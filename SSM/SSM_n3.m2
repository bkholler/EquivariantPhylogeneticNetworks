load("SSM_Networks.m2")

-- setup networks
n = 3
N1 = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 2}}, {4, {1}}})
N2 = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 1}}, {4, {2}}})
N3 = digraph({1,2,3,4,5,6}, {{6, {5, 4, 1}}, {5, {4, 2}}, {4, {3}}})

-- make rings and maps
KK = 
S = KK[ssmStates(n) / (i -> q_i)]
retEdges = {{6, 4}, {5, 4}}
phi1 = ssmFourierNetworkParametrization(retEdges, N1, SourceRing => S, UseStochasticParameters => true);
phi2 = ssmFourierNetworkParametrization(retEdges, N2, SourceRing => S, UseStochasticParameters => true);
phi3 = ssmFourierNetworkParametrization(retEdges, N3, SourceRing => S, UseStochasticParameters => true);

-- load the ideal for N1 which we have already computed
load("SSM_n3_ideal.m2")
-- check that the polynomials in I1 evaluate to zero under phi1
-- this takes about 60 seconds
time all(I1_*, f -> phi1(f) == 0)

-- now we sample points from the variety of N2 and N3
samples1 = sampleImage(10, phi1);
samples2 = sampleImage(10, phi2);
samples3 = sampleImage(10, phi3);

-- the following code shows that 
any(samples1, pt -> sub(I1_0, pt) != 0) -- true which we already know from the above
any(samples2, pt -> sub(I1_0, pt) != 0) -- false, implies M_N1 != M_N2
any(samples3, pt -> sub(I1_0, pt) != 0) -- false, implies M_N1 != M_N3