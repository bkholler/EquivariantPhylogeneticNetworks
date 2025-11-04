load "GroupBasedNetworks.m2"

-- make our map which parametrizes the vanishing ideal of the 4-sunlet
n = 5
M = JCmodel
S = qRing(n, M);
sunletImages = sunletParam(n, M);
R = ring first sunletImages;
psi = map(R, S, sunletImages);


-- make a ring which only keeps the unique variables modulo linear relations
uniqueImPos = apply(unique sunletImages, i -> position(sunletImages, j -> j == i))
KK = QQ
R' = KK[gens R];
S' = KK[(gens S)_uniqueImPos];


-- mod out by the linear reulations
linRels = ideal flatten for i in uniqueImPos list(

	f := sunletImages_i;
	allPos := positions(sunletImages, g -> g == f);
	for j in allPos list if i != j then S_i - S_j else continue
)


-- make our map and run MultigradedImplicitization
-- this should only take a few seconds since we can distinguish with quadratics
phi = map(R', S', apply(sunletImages_uniqueImPos, i -> sub(i, R')));
Hphi = componentsOfKernel(2, phi, UseInterpolation => false);
I = sub(ideal delete(null, flatten values Hphi), S);
J = ideal select(I_*, f -> psi(f) == 0) + linRels;


-- make the map parametrizing the toric variety of the tree
T = leafTree(toList(1..5), {{2, 3}, {1, 2, 3}})
A = phyloToricAMatrix(T, M)
RT = QQ[a_1..a_(numrows(A))]
psiT = map(RT, S, for i to numcols(A)-1 list RT_(flatten entries A_i));


-- get the generators which are not in IT
-- this list is non empty and thus we see IN \not \subset IT as needed
distinguishGens = select(J_*, f -> psiT(f) != 0);
tally(distinguishGens / degree)
f = distinguishGens_0
psiT(f)