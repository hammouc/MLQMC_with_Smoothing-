gen_vec = 4096;
shifts = 1;
M = gen_vec*shifts;
dimension = 3;

%A = RQMC_points(gen_vec, shifts, dimension);
P = sobolset(dimension);
P = scramble(P,'MatousekAffineOwen');
A = transpose(net(P, M));

scatter(A(1,:), A(3,:));