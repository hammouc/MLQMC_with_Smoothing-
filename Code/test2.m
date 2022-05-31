M = 10000;

%Q = randn(M,1);
P = sobolset(1);
P = scramble(P,'MatousekAffineOwen');
Q = norminv(net(P, M));
%% 

estimator = zeros(M,1);
variance = zeros(M,1);
bound = zeros(M,1);
for samples = 1:M
    if Q(samples) > 0
        Q(samples) = 1;
    else
        Q(samples) = 0;
    end
    estimator(samples) = mean(Q(1:samples));
    variance(samples) = std(estimator(1:samples));
    bound(samples) = 1/(samples^0.5);
end
%% 

loglog(1:M, abs(estimator), 1:M, variance, 1:M, bound);