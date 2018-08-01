function y = likelihood_b1_1(X,para)

[N,K] = size(X);
likelihood = 0;

for n = 1:N
    Beta = max(betapdf(X(n,1),para(1),para(2)),0.001);
    Gamma = max(gampdf(X(n,2),para(3),para(4)),0.001);
    Gaussian = max(normpdf(X(n,3),para(5),para(6)),0.001);
    likelihood = likelihood + log(Beta) + log(Gamma) + log(Gaussian);
end

y = -likelihood;

end
