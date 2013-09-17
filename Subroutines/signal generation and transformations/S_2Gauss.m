function [S] = S_2Gauss(size,rho,m_1,m_2,var_1,var_2)
% Creates an homogeneous random Bi Gaussian vector with rho as fraction
% of big components, var_2 as expectation for this gaussian and
% var_2 for it's variance. var_1 is the variance of the small
% components of the produced signal and m_1 its average.
% size is the number of components of the original signal

num_non_zero=floor(rho*size);
num_zero=size-num_non_zero;

SS=([randn(1,num_zero).*sqrt(var_1)+m_1,randn(1,num_non_zero).*sqrt(var_2)+m_2]);
S(randperm(size))=SS;
end




