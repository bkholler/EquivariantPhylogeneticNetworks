load "GroupBasedNetworks.m2"

-- make our map which parametrizes the vanishing ideal of the 4-sunlet
n = 4
M = JCmodel
S = qRing(n, M);
sunletImages = sunletParam(n, M);
R = ring first sunletImages;
psi = map(R, S, sunletImages);

-- make a ring which only keeps the unique variables modulo linear relations
uniqueImPos = apply(unique sunletImages, i -> position(sunletImages, j -> j == i))
--KK = ZZ/nextPrime(32003)
KK = QQ
R' = KK[gens R];
S' = KK[(gens S)_uniqueImPos];

linRels = ideal flatten for i in uniqueImPos list(

	f := sunletImages_i;
	allPos := positions(sunletImages, g -> g == f);
	for j in allPos list if i != j then S_i - S_j else continue
)

-- make our map and run MultigradedImplicitization
phi = map(R', S', apply(sunletImages_uniqueImPos, i -> sub(i, R')));
Hphi = componentsOfKernel(5, phi, UseInterpolation => false);
I = sub(ideal delete(null, flatten values Hphi), S);
J = ideal select(I_*, f -> psi(f) == 0) + linRels;

-- check if our resulting ideal is prime and the right dimension
jac = jacobian matrix psi;
kerDim = max for i from 1 to 100 list rank sub(jac, apply(gens ring jac, i -> i => random(KK)))
rank jacobian matrix phi
isPrime J


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
	all(I_*, f -> psi(f) == 0)
)
