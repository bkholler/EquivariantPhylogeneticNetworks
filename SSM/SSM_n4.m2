needs "SSM_Networks.m2"
needs "../MatroidDistinguish.m2"
needsPackage "MultigradedImplicitization"

-- setup sunlet
n = 4
KK = ZZ/32003;
S = KK[ssmStates(n) / (i -> q_i)]
S4 = digraph(toList(1..8), {{5,1}, {6,2}, {7,3}, {8,4}, {6,5}, {8,5}, {7,8}, {7,6}});
retEdges = {{6, 5}, {8, 5}}

-- compute the map for the 4-sunlet
phi = ssmFourierNetworkParametrization(retEdges, S4, SourceRing => S, UseStochasticParameters => true);
Jphi = jacobian matrix phi;


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


idHash = new MutableHashTable;

-- We know that J \not \subset I so it suffices to check that I \not \subset J
-- we do this by finding a set where Jphi drops rank but Jpsi does not
for i from 0 to 5 list(

    N := networks_i;
    rets := retEdgeSets_i;
    psi := ssmFourierNetworkParametrization(retEdgeSets_0, networks_0, SourceRing => S, UseStochasticParameters => true);
    Jpsi = jacobian matrix psi;
    idHash#N = matroidDistinguish(1000000, 42, Jpsi, Jphi)
)