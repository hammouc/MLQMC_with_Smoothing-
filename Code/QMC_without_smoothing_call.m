d=1;
N=2^3;

gen_vec = 10000;
shifts = 10;
M = gen_vec*shifts;

sigma = 0.4;
T = 1;
K = 100;
X_0 = 100;

%Q = 1/sqrt(N-1).*norminv(RQMC_points(gen_vec, shifts, N-1));
%Q = 1/(N-1)^(0.5).*randn(N-1,M);
P = sobolset(N-1);
P = scramble(P,'MatousekAffineOwen');
Q = 1/sqrt(N-1).*transpose(norminv(net(P, M)));

samplepayoffs = zeros(M,1);
evs = zeros(M,1);

for samples = 1:M
    X = X_0;
    for n = 1:N-1
        X = X + sigma*X*Q(n, samples);
    end
    if X >= K
        samplepayoffs(samples) = X-K;
    end
    evs(samples) = samplepayoffs(samples);
end

refsol = sum(evs)/M

weakerror = zeros(M,1);
rmse = zeros(M,1);
for samples = 1:M
    rmse(samples) = mean(abs(evs(1:samples)-refsol).^2);
    weakerror(samples) = abs(mean(evs(1:samples))-refsol);
end

variance = var(samplepayoffs);
bound = zeros(shifts*gen_vec,1);
boundrate = 0.5;
for j=1:shifts
    for i = 1:gen_vec
        bound(i+(j-1)*gen_vec) = 1.96*sqrt(variance)/(((j-1)*gen_vec + i)^boundrate);
    end
end



loglog(1:M, weakerr, 'blue', 1:M, bound, 'red');
title("QMC for the call option and " + M + " samples (without smoothing)", 'Interpreter', 'latex');
xlabel("M", 'Interpreter', 'latex');
ylabel("Error (Variance: " + variance + ")", 'Interpreter','latex');
legend("Exact error", "Error fit of order $M^{-" + boundrate + "}$", 'Interpreter', 'latex');
%saveas(gcf,'../Slides/Figure/QMC_Call_without_Smoothing.svg');
