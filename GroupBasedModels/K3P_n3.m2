load "GroupBasedNetworks.m2"
needsPackage "Permutations"


-- make our map which parametrizes the vanishing ideal of the 4-sunlet
n = 3
M = K3Pmodel
S = qRing(n, M);
sunletImages = sunletParam(n, M);
R = ring first sunletImages;
psi = map(R, S, sunletImages);


-- make our map and run MultigradedImplicitization
-- this should only take a few seconds since we can distinguish with quadratics
H = componentsOfKernel(4, psi, UseInterpolation => false);
I = sub(ideal delete(null, flatten values H), S);


-- obtain the ideals for the other 2 sunlets by permuting coordinates
(I1, I2, I3) = toSequence for sigma in {{0, 1, 2}, {1, 0, 2}, {2, 1, 0}} list(

    tau := (sigma / (i -> i + 1)) // permutation // inverse // toList / (i -> i -1);
    permuteVars :=  apply(gens S, j -> j => q_(toSequence (toList last baseName j)_tau));
    J := sub(I, permuteVars);
    J
);


-- check that they are all equal
I1 == I2
I1 == I3