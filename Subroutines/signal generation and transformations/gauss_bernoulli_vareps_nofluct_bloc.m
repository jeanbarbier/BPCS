function [S] = gauss_bernoulli_vareps_nofluct_bloc(rho,exp_gauss,sqrt_var_gauss,vareps,L,N1,N2)
% Create an homogeneous random gauss-bernoulli vector with rho as fraction
% of non zero components and exp_gauss as expectation for the gaussian and
% var_gauss for it's variance.
% size : number of components of the original signal
% rho : fraction of non zero components of the original signal

S=[];
for numL=1:L
    if numL==1
        size=N1;
    else
        size=N2;
    end
    num_non_zero=floor(rho*size);
    num_zero=size-num_non_zero;

    SS=([randn(1,num_zero)*sqrt(vareps),randn(1,num_non_zero).*sqrt_var_gauss+exp_gauss]);
    SS_mixed=(intrlv(SS,randperm(size)));
    
    S=[S SS_mixed];

end




