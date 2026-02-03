load "GroupBasedNetworks.m2"

-- make our map which parametrizes the vanishing ideal of the 4-sunlet
n = 4
M = K3Pmodel
S = qRing(n, M);
sunletImages = sunletParam(n, M);
R = ring first sunletImages;
psi = map(R, S, sunletImages);


-- make our map and run MultigradedImplicitization
-- uncomment the following code to compute J as we originally did which takes about 4 hours
-- H = componentsOfKernel(4, psi, UseInterpolation => false);
-- J = sub(ideal delete(null, flatten values H), S);

-- this load the generators of ker(psi) = I_{S_5} which we already computed
use S;
load("K3P_n4_d4_sunlet_gens.m2");


-- we can quickly check that they all do lie in the kernel
all(J_*, f -> psi(f) == 0)


-- all 6 3-sunlet networks with a cherry plus their corresponding reticulation edges
networks = {
digraph({1,2,3,4,5,6,7,8}, {{7, 8}, {7, 6}, {8, 6}, {6, 5}, {5, 1}, {5, 2}, {7, 3}, {8, 4}}),
digraph({1,2,3,4,5,6,7,8}, {{8, 7}, {6, 7}, {8, 6}, {6, 5}, {5, 1}, {5, 2}, {7, 3}, {8, 4}}),
digraph({1,2,3,4,5,6,7,8}, {{7, 8}, {7, 6}, {6, 8}, {6, 5}, {5, 1}, {5, 2}, {7, 3}, {8, 4}}),
digraph({1,2,3,4,5,6,7,8}, {{7, 8}, {7, 6}, {8, 6}, {6, 5}, {5, 3}, {5, 4}, {7, 1}, {8, 2}}),
digraph({1,2,3,4,5,6,7,8}, {{8, 7}, {6, 7}, {8, 6}, {6, 5}, {5, 3}, {5, 4}, {7, 1}, {8, 2}}),
digraph({1,2,3,4,5,6,7,8}, {{7, 8}, {7, 6}, {6, 8}, {6, 5}, {5, 3}, {5, 4}, {7, 1}, {8, 2}}),
}


retEdgeSets = {
{{7, 6}, {8, 6}},
{{6, 7}, {8, 7}},
{{7, 8}, {6, 8}},
{{7, 6}, {8, 6}},
{{6, 7}, {8, 7}},
{{7, 8}, {6, 8}}
}


-- We know that J \not \subset I so it suffices to check that I \not \subset J
for i from 0 to 5 list(

	N := networks_i;
	rets := retEdgeSets_i;
	cycleImages := oneCycleNetworkParam(rets, N, M);
	psi := map(ring first cycleImages, S, cycleImages);
	any(J_*, f -> psi(f) != 0)
)