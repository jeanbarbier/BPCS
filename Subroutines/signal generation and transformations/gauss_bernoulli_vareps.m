function [S] = gauss_bernoulli_vareps(size,rho,exp_gauss,sqrt_var_gauss,vareps)
% Create an homogeneous random gauss-bernoulli vector with rho as fraction
% of non zero components and exp_gauss as expectation for the gaussian and
% var_gauss for it's variance.
% size : number of components of the original signal
% rho : fraction of non zero components of the original signal


    num_non_zero=floor(rho*size);
    num_zero=size-num_non_zero;

    SS=([randn(1,num_zero)*sqrt(vareps),randn(1,num_non_zero).*sqrt_var_gauss+exp_gauss]);
    S=(intrlv(SS,randperm(size)));




