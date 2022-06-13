d=1;
N=2^3;

gen_vec = 10000;
shifts = 8;
M = gen_vec*shifts;

sigma = 0.4;
T = 1;
K = 100;
X_0 = 100;

Q = 1/sqrt(N-1).*norminv(RQMC_points(gen_vec, 200*shifts, N-1));
%Q = 1/(N-1)^(0.5).*randn(N-1,200*M);
%P = sobolset(N-1);
%P = scramble(P,'MatousekAffineOwen');
%Q = 1/sqrt(N-1).*transpose(norminv(net(P, 200*M)));

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
for samples = 1:(200*M)
    X = X_0;
    for n = 1:N-1
        X = X + sigma*X*Q(n, samples);
    end
    if X >= K
        avg = avg + X-K;
    end
end

refsol = avg/(200*M)

weakerror = zeros(M,1);
variances = zeros(M,1);
variance = zeros(M,1);
rmse = zeros(M,1);
boundrate = 0.7;
for samples = 1:M
    variances(samples) = 1/(samples^boundrate);
    weakerror(samples) = abs(mean(evs(1:samples))-refsol);
end

for samples = 1:M
    variance(samples) = std(variances(1:samples));
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
                samplepayoffs(samples) = X-K;
            end
        end
        shiftvar(s) = mean(samplepayoffs);
    end
    gen_shift_var(gen_number) = std(shiftvar);
end




%% 

%variance = var(samplepayoffs);
bound = zeros(M,1);
boundrate = 0.75;
for j=1:M
    bound(j) = 110*1/(j^boundrate);
end

%% 

loglog(1:M, weakerror, gen_vec, 1.96.*gen_shift_var, 1:M, bound);
title("QMC for the call option and " + M + " samples (without smoothing)", 'Interpreter', 'latex');
xlabel("M", 'Interpreter', 'latex');
ylabel("Error (Variance: " + variance(M) + ")", 'Interpreter','latex');
legend("Exact error", "Variance of the estimator (with fitting rate $M^{-" + boundrate + "}$)", 'Interpreter', 'latex');
%saveas(gcf,'../Slides/Figure/QMC_Call_without_Smoothing.svg');
