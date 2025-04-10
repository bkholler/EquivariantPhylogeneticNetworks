load("../../PhylogeneticNetworks.m2")
load("../../MatroidDistinguish.m2")
load("../GMM_Networks.m2")
needsPackage "NAGtypes"

-- create the ring for 2-state GMM on a 4-leaf network
n = 4
k = 3
S = pRing(k, n)

-- setup grading
A = transpose matrix for x in gens(S) list(

	l = last baseName x; 

	flatten for l' in l list apply(k, i -> if i == l' then 1 else 0)
);

B = submatrix'(A, toList(0..k-1),)
S = newRing(S, Degrees => B);


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
I1 = minors(k, flat12');
I2 = minors(k, flat14');

-- look at all possible products of these minors
prods = flatten for f in I1_* list for g in I2_* list f*g;
polysByDeg = hashTable(join, apply(prods, i -> {degree(i), {i}}));


-- filter the lattice points for supports which are actually possible circuits
goodLats = for deg in keys(polysByDeg) list(

    supp := support matrix {polysByDeg#deg};

    if rank(jac_(supp / index)) < #supp then deg else continue
);

tally apply(goodLats, i -> #(polysByDeg#i))
remainLats = new MutableHashTable from apply(goodLats, i -> {i, {}})



-- filter the result for polynomials which actually lie in the vanishing ideal
load("gmm3_point_sample.m2");
KK = ZZ/32003;
pts = V / (v -> sub(matrix {values hashTable v}, KK));

for deg in keys(remainLats) list(

    componentSize = #polysByDeg#deg;
    B = matrix{polysByDeg#deg};

    M = matrix for pt in take(pts, componentSize) list flatten entries sub(evaluate(B, pt), KK);
    curRank = rank(M);

    if curRank < componentSize then(
        H#deg = M;
        print(concatenate("!!Found Poly in!!: ", toString(deg)));
        remove(remainLats, deg);
    );

    if curRank == componentSize then(
        print(concatenate("no polys found in: ", toString(deg)));
        remove(remainLats, deg);
    )
)
