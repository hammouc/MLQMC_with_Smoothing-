function A = RQMC_points(genvec_len,M_qmc,N)
%   Produces a N x genvec_len*M_qmc matrix with QMC numbers / N to be the
%   dimension for example
    generators = i4_sobol_generate(N,genvec_len,0); % NxM Matrix
    shifts = rand(N,M_qmc);
    A = zeros(N,M_qmc*genvec_len);
    for s = 1:M_qmc
        shifted_points = mod(shifts(:,s) + generators, 1);
        A(:,((s-1)*genvec_len+1):(s*genvec_len)) = shifted_points;
    end
end

