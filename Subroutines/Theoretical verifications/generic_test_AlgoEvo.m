disp('Test')

% signal and test properties
% choice of the test
choice = 'L1';
% proportion of measurements
measure_rate = alpha;
% size of the signal..
size_S = 20000;
% ..and it's density
rho = 1;
% bloc matrix or not?
bloc = 0;

% creation of the signal
switch choice
    case 'SparseGauss'
        mean = 0; var = 1; S = S_SparseGauss(size_S,rho,mean,var);
    case 'SparseExponential'
        exposent = 3; S = S_SparseExponential(size_S,rho,exposent);
    case 'SparseBinary'
        S = S_SparseBinary(size_S,rho);
    case '2Gauss'
        m_2 = 0; var_1 = 1e-7; var_2 = 2; S = S_2Gauss(size_S,rho,m_2,var_2,var_1);
    case 'SparseConstant'
        down = -1; up = 1;  S = S_SparseConstant(size_S,rho,down,up);
    case 'L1'
        %         min = 0; max = 1; S = S_SparseConstant(size_S,rho,min,max);
        %         mean = 0; var = 1; S = S_SparseGauss(size_S,rho,mean,var);
        %         S = S_SparseBinary(size,rho);
        S = (-log(rand(size_S,1)).*(2.*(rand(size_S,1)>0.5)-1))';
        %         beta = 1; S = 1 ./ beta .* log(rand(1,size_S) .* beta ./ 2) .* sign(2 .* (rand(1,size_S) ) - 1);
    case 'Laplace'
        %         min = 0; max = 1; S = S_SparseConstant(size_S,rho,min,max);
        %         mean = 0; var = 1; S = S_SparseGauss(size_S,rho,mean,var);
        %         S = S_SparseBinary(size_S,rho);
        %        beta = 1; S = 1 ./ beta .* log(rand(1,size_S) .* beta ./ 2) .* sign(2 .* (rand(1,size_S) ) - 1);
        S = (-log(rand(size_S,1)).*(2.*(rand(size_S,1)>0.5)-1))';
    otherwise
        disp('Error in choice')
end


% creation of the measurement matrix
if (bloc == 1)
    L = 30; J1 = 1 ./ size_S; J2 = 10 ./ size_S; alpha1 = 0.5; alpha2 = 0.4;
    [G,Jmat,Mvec,Nvec] = Create_MAT1D(size_S,L,J1,J2,alpha1,alpha2);
    %     J = 1;
    %     [F,Jmat,Mvec,Nvec] = Create_MAT1D_gauss(size_S,L,J,alpha1,alpha2);
    G = G ./ sqrt(size_S);
    Jmat = Jmat ./ size_S;
    Jfull = convert2fullMN(Jmat,Mvec,Nvec);
else
    G = randn(floor(measure_rate .* size_S),size_S) ./ sqrt(size_S);
end

% average of the matrix
mean_G = 0.;
G = G + mean_G ./ sqrt(size_S);
% measurement noise
var_noise = 0;
% the measurement
Y = G * S' + sqrt(var_noise) .* randn(min(size(G) ),1);

% algorithm properties
My = CSBP_Solver_Opt();
My.nb_iter = 200;
My.save_speed = 1;
My.save_memory = 1 - My.save_speed;
My.print = 1;
My.conv = 10^(-15);
My.learn = 0;
My.learn_d = 0;
My.signal_rho = -1; % (-1 -> measure_rate ./ 10)
My.var_noise = 0;
My.signal = S;
My.dump_learn = 0.;
My.dump_mes = 0.5;
My.prior = choice;
My.option_noise = 0;
My.remove_mean = 0;
My.method = 'AMP';
% My.varG = Jfull;
% My.Nvec_bm = Nvec;
% My.Mvec_bm = Mvec;
% SparseGauss
My.m_gauss = 0;
My.var_gauss = 1;
% SparseExponential
My.expo = 1;
% SparseConstant
My.c_down = -1;
My.c_up = 1;
% 2Gauss
My.m_2_gauss = 0;
My.var_1_gauss = 1e-10;
My.var_2_gauss = 1;
% L1
My.min = -15;
My.max = 15;
% Laplace
My.beta = 1;

% reconstruction
tic;
[results,n_and_e] = CSBP_Solver(Y,G,My);
toc;