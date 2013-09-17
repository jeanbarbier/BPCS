% EMBP_partial_Opt
% Function to set EMBP_Solver_Opt to their default values
%
%    Details of the option:
%   .method             method of resolution of the inference problem (see
%                       the article by Krzakala and al.). It can be 'AMP', 'AMPh' for statistically 
%                       homogeneous measurements matricesor, 'BP' or AMPHb for bloc matrices. ['AMP']
%   .save_memory        option that does not allow to create new objects of the
%                       size comparable with the measurement matrix during
%                       the procedure, such that the memory required by the
%                       routine is lowered. [1]
%   .save_speed         option that allows pre-computation of new objects of the
%                       size comparable with the measurement matrix, such that the memory required by the
%                       routine is multiplied by a factor 3 to 4 times the one required to store the measurement matrix
%                       but the algorithm is speeded-up. It is usefull for not too big matrices or high memory capacities. [0]
%                       WARNING: save_memory and save_speed are not
%                       compatible: One must be equal to 1, the other to 0.
%   .nb_iter            max number of iterations [default : 1000].
%   .print              print results every opt.print iterations (0 -> never) [10].
%   .conv               convergence criterion [1e-8].
%   .learn              learning parameters or not [0].
%   .learn_d            learning dumping or not [0].
%   .signal_rho         first estimate of the density of the signal [M/10N].
%   .var_noise          first estimate of the variance vector of the noise [1e-10].
%   .signal             a solution to compare to while running.
%   .dump_learn         dumping coefficient of the learning [0.].
%   .dump_mes           dumping of the messages [0.5].
%   .option_noise       Value is 1/~=1 [0, i.e not activated]. If activated,
%                       it deals with matrix uncertainty.
%   .remove_mean        [0] When activited (1 or 2 instead of 0), this
%                       allows the algorithm to deal with non zero
%                       mean matrices that can have a different mean on every column (1).
%                       With value 2 it assumes the same average value for each of them.
%                       If provided (see bellow), it uses the one given by the user.
%   .Gmean              see .remove_mean
%   .Ymean              see .remove_mean
%   .Nvec               [1] This can be used to specify a structure
%                       to the signal, and change the printing output.
%   .Mvec               [1] This can be used to specify a structure
%                       to the solution vector Y (this is used in
%                       Active option.noise for instance).
%   .varG               Used if AMPhb is activated (for bloc matrices). It
%                       can be given in two different forms: Or it is a full sparse [M N] matrix where the element varG(i,j) is the variance of
%                       the bloc to which G(i,j) belong to, or it a small [L C]
%                       matrix where the element varG(i,j) is the variance
%                       of the bloc (i,j). In the second case, the user
%                       must also give the opt.Mvec_bm and opt.Nve_bm. 
%   .Nvec_bm            Used when AMPhb is activated. It is a vector containing the number of columns in each bloc. 
%   .Mvec_bm            Used when AMPhb is activated. It is a vector containing the number of lines in each bloc. 
%   .prior              prior on the data ['SparseGauss'].
%                       One can use 'SparseGauss', 'SparseBinary',
%                       '2Gauss' for approximate-sparsity with 2
%                       gaussians, 'SparseExponential', 'SparseConstant', 'L1' and 'Laplace'.
%   .m_gauss            Mean of the gaussian for the prior 'SparseGauss' [0].
%   .var_gauss          Variance .. [1.]
%   .expo               Exposent for the prior 'SparseExponential' [1];
%   .c_down             Lower bound of the signal for the prior 'SparseConstant' [0]
%   .c_up               Upper bound .. [1]
%   .m_2_gauss          Mean of the second gaussian in the prior, so the one of big components of density rho [0].
%   .var_2_gauss        Variance .. [1]
%   .var_1_gauss        Variance of the small gaussian of mean 0 that models the small components of density (1 - rho) [1e-5].
%   .min                Minimum bound for the L1 reconstruction [-3].
%   .max                Max bound for the L1 reconstruction [3].
%   .beta               inverse temperature for the Laplace prior [1].

function opt = CSBP_Solver_Opt()
opt.nb_iter = 1000;
opt.print = 10;
opt.conv = 10^(-8);
opt.learn = 0;
opt.learn_d = 0;
opt.signal_rho = -1;
opt.var_noise = 10^(-10);
opt.signal = [];
opt.dump_learn = 0.;
opt.dump_mes = 0.5;
opt.prior = 'SparseGauss';
opt.save_memory = 1;
opt.save_speed = 0;
opt.option_noise = 0;
opt.remove_mean = 0;
opt.Gmean = [];
opt.Ymean = [];
opt.Nvec = [1];
opt.Mvec = [1];
opt.varG = [];
opt.Nvec_bm = [];
opt.Mvec_bm = [];
opt.method = 'AMP';
opt.alphaBig = 0;
opt.MSEbyBlock = 0;
% SparseGauss
opt.m_gauss = 0;
opt.var_gauss = 1;
% SparseExponential
opt.expo = 1;
% SparseConstant
opt.c_down = 0;
opt.c_up = 1;
% 2Gauss
opt.m_1_gauss = 0;
opt.m_2_gauss = 0;
opt.var_1_gauss = 1e-5;
opt.var_2_gauss = 1;
% L1
opt.min = -3;
opt.max = 3;
% Laplace
opt.beta = 1;

end