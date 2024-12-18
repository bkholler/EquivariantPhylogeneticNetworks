N = digraph({1,2,3,4,5,6}, {{6, {5, 4, 3}}, {5, {4, 2}}, {4, {1}}})
k = 2

-- since the jacobian is full rank, the resulting variety is the whole space
phi = gmmNetParam(k, 6, N, {{5, 4}, {6, 4}})
J = jacobian matrix {phi}
rank(J)