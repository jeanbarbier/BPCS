clear all

% signal and test properties
% choice of the test
choice = 'Binary1';
% proportion of measurements
measure_rate = 0.55;
% size of the signal..
size_S = 10000;
% ..and it's density
rho = .1;
% bloc matrix or not?
bloc = 0;

% creation of the signal
switch choice
    case 'SparseGauss'
        average = 0; var = 1; S = S_SparseGauss(size_S,rho,average,var);
    case 'SparseGaussCut'
        average = 0; var = 1; S = S_SparseGauss(size_S,rho,average,var);
    case 'SparseExponential'
        exposent = 3; S = S_SparseExponential(size_S,rho,exposent);
    case 'SparseBinary'
        S = S_SparseBinary(size_S,rho);
    case '2Gauss'
        m_2 = 0; var_1 = 1e-6; var_2 = 1; S = S_2Gauss(size_S,rho,m_2,var_2,var_1);
    case 'SparseConstant'
        down = -1; up = 1;  S = S_SparseConstant(size_S,rho,down,up);
    case 'L1'
        %         min = 0; max = 1; S = S_SparseConstant(size_S,rho,min,max);
        %         mean = 0; var = 1; S = S_SparseGauss(size_S,rho,mean,var);
        %         S = S_SparseBinary(size,rho);
        S = (-log(rand(size_S,1)).*(2.*(rand(size_S,1)>0.5)-1))';
    case 'Laplace'
        %         min = 0; max = 1; S = S_SparseConstant(size_S,rho,min,max);
        %         mean = 0; var = 1; S = S_SparseGauss(size_S,rho,mean,var);
        %         S = S_SparseBinary(size_S,rho);
        S = (-log(rand(size_S,1)).*(2.*(rand(size_S,1)>0.5)-1))';
    case 'Binary1'
        S = S_SparseBinary(size_S,0.5);
        S = S - (S == 0);
        rho = 1;
    otherwise
        disp('Error in choice')
end


% creation of the measurement matrix
if (bloc == 1)
    %     M = ceil(size_S .* measure_rate);
    %     L = 4; J1 = 1 ./ size_S; J2 = 0. ./ size_S;
    %     [G,Jmat,Mvec,Nvec] = Create_MAT1D_RC(size_S,M,L,J1,J2,alpha2);
    %         J = 1;
    %         [G,Jmat,Mvec,Nvec] = Create_MAT1D_gauss(size_S,L,J,alpha1,alpha2);
    %     G = G ./ sqrt(size_S);
    %     Jmat = Jmat ./ size_S;
    %     Jfull = convert2fullMN(Jmat,Mvec,Nvec);
    L = 10; J = .2; alpha1 = .22; alpha2 = .16; W = 3;
    [G,Jmat,Mvec,Nvec] = Create_MAT1D_gauss_Wplus(size_S,L,J,alpha1,alpha2,W);
    %     [G,Jmat,Mvec,Nvec] = Create_MAT1D_gauss_Fullplus(size_S,L,J,alpha1,alpha2);
else
    G = randn(ceil(measure_rate .* size_S),size_S) ./ sqrt(size_S);
end

% average of the matrix
mean_G = 0.;
G = G + mean_G ./ sqrt(size_S);
% measurement noise
var_noise = 0.01;
% the measurement
Y = G * S' + sqrt(var_noise) .* randn(min(size(G) ),1);

% algorithm properties
My = CSBP_Solver_Opt();
My.nb_iter = 1000;
My.save_speed = 1;
My.save_memory = 1 - My.save_speed;
My.print = 1;
My.conv = 0;
My.learn = 0;
My.learn_d = 0;
My.signal_rho = -1; % (-1 -> measure_rate ./ 10)
if (My.learn == 0); My.signal_rho = rho; end;
My.calta = 0;
My.var_noise = 0;
My.signal = S;
My.m_S = mean(S);
My.m_S2 = mean(S.^2);
My.dump_learn = 0.;
My.dump_mes = 0.5;
My.prior = choice;
My.option_noise = 0;
My.remove_mean = 0;
My.method = 'AMP';
My.alphaBig = 0;
if (bloc ==1)
    %     My.varG = Jfull;
    %     My.Nvec_bm = Nvec;
    %     My.Mvec_bm = Mvec;
end
% SparseGauss
My.m_gauss = 0;
My.var_gauss = 1;
% SparseGaussCut
My.mC_gauss = 0;
My.varC_gauss = 1;
My.Cut = 5;
% SparseExponential
My.expo = 1;
% SparseConstant
My.c_down = -1;
My.c_up = 1;
% 2Gauss
if (strcmp(choice,'2Gauss') )
    My.m_2_gauss = m_2;
    My.var_1_gauss = var_1;
    My.var_2_gauss = var_2;
end
% L1
My.min = -10;
My.max = 10;
% Laplace
My.beta = 1;

% reconstruction
tic;
[results,n_and_e] = CSBP_Solver(Y,G,My);
toc;