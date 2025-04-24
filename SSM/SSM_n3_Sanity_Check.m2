load("SSM_Networks.m2")

-- setup networks
n = 3
retEdges = {{6, 4}, {5, 4}}
N1 = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 2}}, {4, {1}}})
N2 = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 1}}, {4, {2}}})
N3 = digraph({1,2,3,4,5,6}, {{6, {5, 4, 1}}, {5, {4, 2}}, {4, {3}}})

-- make the parametrization in the probability coordinates
P = pRing(4, n)
phi1 = ssmNetworkReParametrization(retEdges, N1, UseStochasticParameters => true);
phi2 = ssmNetworkReParametrization(retEdges, N2, UseStochasticParameters => true);
phi3 = ssmNetworkReParametrization(retEdges, N3, UseStochasticParameters => true);

-- load the ideal for N1 which we have already computed
KK = QQ;
S' = KK[p_(toList(n:{0,0}))..p_(toList(n:{1,1}))];
S = KK[q_(toList(n:{0,0}))..q_(toList(n:{1,1}))];
load("SSM_n3_ideal.m2")
I1 = sub(I1, S);

-- sample points from the variety of N1, N2, and N3
samples1 = sampleImage(10, phi1);
samples2 = sampleImage(10, phi2);
samples3 = sampleImage(10, phi3);

-- transform I1 into probability coordinates
probL = for x in gens(S') list(

    ind := last baseName x;

    j := apply(ind, i -> i_0);
    i := apply(ind, i -> i_1);

    sum for k in toList((n:0)..(n:1)) list(

        sign := (-1)^(sum toList apply(j, k, (jl, kl) -> jl*kl));
        pInd := apply(k, i, (kl, il) -> {kl, il});

        sign*p_pInd
    )
);

(M, H) := (coefficients matrix {probL});
H = sub(H, KK);
invFourierTransform = map(P, S, matrix({gens P})*H);
-- note that in probability coordinates the polynomials both involve all 64 coordinates
-- this means the matroid could not distinguish these ideals
J1 = invFourierTransform(I1);


any(samples1, pt -> sub(J1_0, pt) != 0) -- false, this shows every sample is zero 
any(samples2, pt -> sub(J1_0, pt) != 0) -- true, implies M_N1 != M_N2
any(samples3, pt -> sub(J1_0, pt) != 0) -- true, implies M_N1 != M_N3