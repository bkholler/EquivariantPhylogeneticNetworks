load("../../PhylogeneticNetworks.m2")
load("../GMM_Networks.m2")
needsPackage "MultigradedImplicitization"

-- create the ring and the vanishing ideal of the 2-state GMM on a 4-sunlet with reticulation 1 and leaves labelled cyclically
n = 5
k = 2
S = pRing(k, n)


-- compute some polynomials in the vanishing ideal of the n-cycle by interpolating minors
N = digraph(toList(1..10), {{8, 9}, {8, 7}, {9, 10}, {10, 6}, {7, 6}, {6, 1}, {7, 2}, {8, 3}, {9, 4}, {10, 5}})
phiN = gmmNetworkReParametrization(k, {{10, 6}, {7, 6}}, N);


-- remake the rings over a finite field so we can do interpolation
KK = ZZ/32003
R = KK[gens target phiN]
S' = KK[gens source phiN]
psiN = map(R, S', sub(matrix phiN, R));


-- grade the domain
A = transpose matrix for x in gens(S) list(

	l = last baseName x; 

	flatten for l' in l list apply(k, i -> if i == l' then 1 else 0)
);
B = submatrix'(A, toList(0..k-1),)
S = newRing(S, Degrees => B)


-- make the flattening matrices whose minors we will interpolate
flat12 = flat({1, 2}, {3, 4, 5}, S)
flat15 = flat({1, 5}, {2, 3, 4}, S)
minors12 = (minors(3, flat12))_*;
minors15 = (minors(3, flat15))_*;
degs12 = minors12 / degree;
degs15 = minors15 / degree;


-- refine their products by multidegree
detsByDegs = hashTable(join, flatten for i to #minors12 - 1 list for j to #minors15 - 1 list {degs12_i + degs15_j, {{i, j}}});


-- sample points in the image of the parametrization
pts = sampleImage(100, psiN);
remainDegs = new MutableHashTable;
goodDegs = new MutableHashTable;


-- run through our different sets of minors
for deg in keys(detsByDegs) do(

    if #(detsByDegs#deg) > #pts then(
        remainDegs#deg = detsByDegs#deg;
        continue;
    );

    G := detsByDegs#deg / (l -> minors12_(l_0)*minors15_(l_1));

    evalG := matrix for i from 0 to #G - 1 list for j from 0 to #G - 1 list sub(G_j, pts_i);

    if rank(evalG) < #G then(
        goodDegs#deg = evalG;
        --print("found good deg");
    );
)


-- make our ideal of new generators
I = ideal flatten for k in keys(goodDegs) list(

	M := goodDegs#k;
	K := sub(gens ker M, QQ);
	polys = matrix {detsByDegs#k / (l -> minors12_(l_0)*(minors15)_(l_1))};
	flatten entries(polys*K)
);




-- select the generators we found which are not in the ideal of the tree
T = digraph(toList(1..8), {{8, 6}, {8, 7}, {8, 1}, {6, 2}, {6, 3}, {7, 4}, {7, 5}})
phiT = gmmTreeParametrization(k, T, SourceRing => S);

-- Despite the name, when probabilisticMembershipTest returns false then it is true that f \notin \ker(\phi) symbolically
-- the output is only probabilistic when checking that it does belong to the ideal
distinguishGens = select(I_*, f -> not probabilisticMembershipTest(10, f, phiT));

-- check one polynomial symbolically
f = distinguishGens_0
phiN(sub(f, source phiN)) == 0
