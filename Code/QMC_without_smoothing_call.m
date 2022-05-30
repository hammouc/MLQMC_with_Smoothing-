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
Q = 1/sqrt(N-1).*transpose(norminv(net(P, 15*M)));

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

avg = 0;
for samples = 1:(15*M)
    X = X_0;
    for n = 1:N-1
        X = X + sigma*X*Q(n, samples);
    end
    if X >= K
        avg = avg + X-K;
    end
end

refsol = avg/(15*M);

weakerror = zeros(M,1);
variances = zeros(M,1);
variance = zeros(M,1);
rmse = zeros(M,1);
for samples = 1:M
    variances(samples) = mean(evs(1:samples));
    variance(samples) = std(variances(1:samples));
    weakerror(samples) = abs(mean(evs(1:samples))-refsol);
end
%% 

%variance = var(samplepayoffs);
bound = zeros(shifts*gen_vec,1);
boundrate = 0.5;
for j=1:M
    bound(j) = 814*1/j^boundrate;
end

%% 

loglog(1:M, weakerror, 'blue', 1:M, variance, 'red', 1:M, bound);
title("QMC for the call option and " + M + " samples (without smoothing)", 'Interpreter', 'latex');
xlabel("M", 'Interpreter', 'latex');
ylabel("Error (Variance: " + variance(M) + ")", 'Interpreter','latex');
legend("Exact error", "Variance of the estimator (with fitting rate $M^{-" + boundrate + "}$)", 'Interpreter', 'latex');
%saveas(gcf,'../Slides/Figure/QMC_Call_without_Smoothing.svg');
