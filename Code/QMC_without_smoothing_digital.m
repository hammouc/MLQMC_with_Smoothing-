d=1;
N=2^3;

gen_vec = 10000;
shifts = 1;
M = gen_vec*shifts;

sigma = 0.4;
T = 1;
K = 100;
X_0 = 100;

Q = 1/sqrt(N-1).*norminv(RQMC_points(gen_vec, shifts, N-1));
%Q = 1/(N-1)^(0.5).*randn(N-1,M);
%P = sobolset(N-1);
%P = scramble(P,'MatousekAffineOwen');
%Q = 1/sqrt(N-1).*transpose(norminv(net(P, M)));
samplepayoffs = zeros(M,1);
evs = zeros(M,1);

avg = 0;
rmse = zeros(M,1);
for samples = 1:M
    X = X_0;
    for n = 1:N-1
        X = X + sigma*X*Q(n, samples);
    end
    if X > K
        samplepayoffs(samples) = 1;
    end
end

refsol = sum(samplepayoffs)/M;

weakerror = zeros(M,1);
bound = zeros(M,1);
boundrate = 0.5;

for samples=1:M
    weakerror(samples) = abs(mean(samplepayoffs(1:samples))-refsol);
    bound(samples) = 0.6/(samples^boundrate);
end
%% 
d=1;
N=2^3;

gen_vec = pow2(4:16);
shifts = 8;

sigma = 0.4;
T = 1;
K = 100;
X_0 = 100;

%Q = 1/(N-1)^(0.5).*randn(N-1,200*M);
%P = sobolset(N-1);
%P = scramble(P,'MatousekAffineOwen');
%Q = 1/sqrt(N-1).*transpose(norminv(net(P, 200*M)));

gen_shift_var = zeros(numel(gen_vec),1);

for gen_number = 1:numel(gen_vec)
    Q = 1/sqrt(N-1).*norminv(RQMC_points(gen_vec(gen_number), shifts, N-1));
    shiftvar = zeros(shifts,1);
    for s = 1:shifts
        samplepayoffs = zeros(gen_vec(gen_number),1);
        for samples = 1:gen_vec(gen_number)
            X = X_0;
            for n = 1:N-1
                X = X + sigma*X*Q(n, samples+(s-1)*gen_vec(gen_number));
            end
            if X >= K
                samplepayoffs(samples) = 1;
            end
        end
        shiftvar(s) = mean(samplepayoffs);
    end
    gen_shift_var(gen_number) = std(shiftvar);
end

%% 

loglog(1:M, weakerror, 'blue', gen_vec, 1.96*gen_shift_var, 1:M, bound);
title("QMC for the digital option and " + M + " samples (without smoothing)", 'Interpreter', 'latex');
xlabel("M", 'Interpreter', 'latex');
ylabel("Error (Variance: " + variance(M) + ")", 'Interpreter','latex');
legend("Exact error", "Error fit of order $M^{-" + boundrate + "}$", 'Interpreter', 'latex');
%saveas(gcf,'../Slides/Figure/QMC_Digital_without_Smoothing.svg');
