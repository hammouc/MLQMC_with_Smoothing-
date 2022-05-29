d=1;
N=2^2;

gen_vec = 100000;
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
for samples = 1:M
    X = X_0;
    for n = 1:N-1
        X = X + sigma*X*Q(n, samples);
    end
    if X > K
        samplepayoffs(samples) = 1;
    end
    avg = avg + samplepayoffs(samples);
    evs(samples) = avg/samples;
end

refsol = evs(M)
variance = var(samplepayoffs);
bound = zeros(M,1);
boundrate = 0.7;
for i=1:M
    bound(i) = 1.96*sqrt(variance)/(i^boundrate);
end

loglog(1:M, abs(evs-refsol), 'blue', 1:M, bound, 'red');
title("QMC for the digital option and " + M + " samples (without smoothing)", 'Interpreter', 'latex');
xlabel("M", 'Interpreter', 'latex');
ylabel("Error (Variance: " + variance + ")", 'Interpreter','latex');
legend("Exact error", "Error fit of order $M^{-" + boundrate + "}$", 'Interpreter', 'latex');
%saveas(gcf,'../Slides/Figure/QMC_Digital_without_Smoothing.svg');
