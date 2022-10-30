Q = eye(2)

A = [0, 1; -0.5, -0.5]

C = eye(4) - kron(A.', A.')

vec_Q = reshape(Q,1,[]).'

vec_P = inv(C)*vec_Q

P = reshape(vec_P, 2,2)