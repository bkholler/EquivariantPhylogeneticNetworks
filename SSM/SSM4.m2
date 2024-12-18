load "SSM_Networks.m2"

n = 4;
S4 = QQ[ssmStates(n) / (i -> q_i)];
flat12 = flat(n, {1,2}, {3,4}, S4);
flat14 = flat(n, {1,4}, {2,3}, S4);

inds = {{({0, 1, 2},{0, 1, 2}), ({1, 2, 3},{1, 2, 3})}, {({0, 1, 2},{0, 1, 3}), ({0, 2, 3},{1, 2, 3})}, {({0, 1, 2},{0, 2, 3}),
      ({1, 2, 3},{0, 2, 3})}, {({0, 1, 2},{1, 2, 3}), ({0, 2, 3},{0, 2, 3})}, {({0, 2, 3},{0, 1, 2}), ({0, 1, 3},{1, 2, 3})},
      {({0, 2, 3},{0, 1, 3}), ({0, 1, 2},{1, 2, 3})}, {({0, 2, 3},{0, 2, 3}), ({0, 1, 3},{0, 2, 3})}, {({0, 2, 3},{1, 2, 3}), ({0,
      1, 2},{0, 2, 3})}};

coeffs = {1, -1, 1, -1, -1, 1, -1, 1};


f = fillTemplate(inds, coeffs, flat12, flat14, S4);

maps = new MutableHashTable;

for N in gen4Cycles() do(

  images = ssmNetParam(N, {{6,5}, {8,5}});
  curR = ring images_0;

  maps#N = map(curR, S4, images)
  );


N = (gen4Cycles())_0
phi = maps#N;


idVals = flatten flatten for e in {{6, 2}, {7, 3}, {8, 4}} list(
  for j from 0 to 1 list(
    for i in (0,0)..(1,1) list(

      if sum(toList i) % 2 == 1 then a_(e, j, i_0, i_1) => 0 else a_(e, j, i_0, i_1) => 1
    )
  )
)


simplePhi = map(target phi, S4, sub(matrix phi, idVals));

imf = simplePhi(f);


R = QQ[e_({0,0},{0,0})..e_({1,1},{1,1})]
H = new MutableHashTable;






for N in keys(maps) list {certifyPoly(10, f, maps#N), certifyPoly(10, g, maps#N)}



