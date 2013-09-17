function [S] = S_SparseExponential(size,rho,alpha)
% Creates an exponential-Bernoulli distributed signal S of density rho of components with law p(x_i) ~ exp(-alpha * s_i) * I(s_i > 0) 
% where I() is the indicator function.

SS = [zeros(1,size - ceil(rho .* size)),1 ./ alpha .* exprnd(1 ./ alpha,1,ceil(rho .* size))];
S = (intrlv(SS,randperm(size)));

end

