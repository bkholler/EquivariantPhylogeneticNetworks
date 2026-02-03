load "GroupBasedNetworks.m2"

-- make our map which parametrizes the vanishing ideal of the 4-sunlet
n = 5
M = K3Pmodel
S = qRing(n, M);
sunletImages = sunletParam(n, M);
R = ring first sunletImages;
psi = map(R, S, sunletImages);


-- make our map and run MultigradedImplicitization
-- this should only take a few seconds since we can distinguish with quadratics
H = componentsOfKernel(2, psi, UseInterpolation => false);
I = sub(ideal delete(null, flatten values H), S);


-- make the map parametrizing the toric variety of the tree
T = leafTree(toList(1..5), {{2, 3}, {1, 2, 3}})
A = phyloToricAMatrix(T, M)
RT = QQ[a_1..a_(numrows(A))]
psiT = map(RT, S, for i to numcols(A)-1 list RT_(flatten entries A_i));


-- get the generators which are not in IT
-- this list is non empty and thus we see IN \not \subset IT as needed
distinguishGens = select(I_*, f -> psiT(f) != 0);
tally(distinguishGens / degree)
f = distinguishGens_0
psiT(f)