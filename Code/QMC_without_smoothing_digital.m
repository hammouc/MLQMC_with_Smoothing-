d=1;
N=2^3;

gen_vec = 10000;
shifts = 10;
M = gen_vec*shifts;

sigma = 0.4;
T = 1;
K = 100;
X_0 = 100;

Q = 1/sqrt(N-1).*norminv(RQMC_points(gen_vec, shifts, N-1));
%Q = 1/(N-1)^(0.5).*randn(N-1,M);
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
%% 

refsol = sum(samplepayoffs)/M
variance = zeros(M,1);
bound = zeros(M,1);
weakerror = zeros(M,1);
variances = zeros(M,1);
boundrate = 1;
for samples=1:M
    variances(samples) = mean(samplepayoffs(1:samples));
    variance(samples) = var(variances(1:samples));
    weakerror(samples) = abs(mean(samplepayoffs(1:samples))-refsol);
    bound(samples) = 1.96*variance(100)/(samples^boundrate);
end
%% 

loglog(1:M, weakerror, 'blue', 1:M, variance, 'red', 1:M, bound);
title("QMC for the digital option and " + M + " samples (without smoothing)", 'Interpreter', 'latex');
xlabel("M", 'Interpreter', 'latex');
ylabel("Error (Variance: " + variance(M) + ")", 'Interpreter','latex');
legend("Exact error", "Error fit of order $M^{-" + boundrate + "}$", 'Interpreter', 'latex');
%saveas(gcf,'../Slides/Figure/QMC_Digital_without_Smoothing.svg');
