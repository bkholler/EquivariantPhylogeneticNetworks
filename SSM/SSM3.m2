load "SSM_Networks.m2"
needsPackage "Visualize"
needsPackage "MultigradedImplicitization"
needsPackage "NumericalImplicitization"

KK = ZZ/nextPrime(10000000)
FF = QQ

-- If we contract the leaf edge on the reticulation, do we get the same dimension?
-- Can we find the degree of this hypersurface and then compute it?
N = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 2}}, {4, {1}}})
phi = ssmNetParam(N, {{5, 4}, {6, 4}});
R = FF[gens ring phi_0]



n = 3
S = FF[ssmStates(n) / (i -> q_i)];
psi = map(R, S, apply(phi, i -> sub(i, R)));


A = maxGrading psi;
dom = newRing(S, Degrees => A);
B = basis(4, S);
lats = apply(flatten entries sub(B, dom), i -> degree i) // unique;

basisHash = new MutableHashTable;
scan(lats, deg -> basisHash#deg = basis(deg, dom));
gensHash = new MutableHashTable;


FF = ZZ/nextPrime(1000000)
jac = jacobian matrix psi;
jac = sub(jac, apply(gens target psi, t -> t => random(FF)));



count = 0;
skips = 0;

goodLats = new MutableHashTable;


-- This loop finds the lattice points which might yield polynomials
for deg in lats do (
      
  C = findSupportIndices(support sub(basisHash#deg, source psi), psi);

  if count% 100 == 0 then print({count, skips}); 

  if (numcols(basisHash#deg) == 1) then(

    gensHash#deg = {};
    count = count+1;
    skips = skips + 1;
    continue;
    );

  if rank(jac_C) == #C then(

      gensHash#deg = {};
      count = count+1;
      skips = skips + 1;
      continue;
      );

  goodLats#deg = basisHash#deg;

  count = count + 1;

  );


deg = (keys goodLats)_0
B = flatten entries basisHash#deg;

samplePts = for i from 0 to (#B-1) list(

	params = apply(gens target psi, i -> i => random(FF));
	pt = apply(gens source psi, j -> sub(j, dom) => sub(psi(j), params))
);


M = matrix for pt in samplePts list for m in B list sub(m, pt)


f = ((basisHash#deg)*(sub(gens ker M, dom)))_(0,0);

certifyPoly(10000, f, psi)



N' = digraph({1,2,3,4,5,6}, {{6, {5, 4, 2}}, {5, {4, 3}}, {4, {1}}})
imagesPsi = ssmNetParam(N', {{5, 4}, {6, 4}});
psi = map(ring imagesPsi_0, S, imagesPsi);
certifyPoly(10000, f, psi)