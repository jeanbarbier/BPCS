addpath(genpath('./BPCS'));

%% Parameters
% size of the phantom
N = 32.^2;
% seeded or homogeneous matrix (if seeded = 0)?
seeded =1;
% variance of the small additive gaussian noise z and density of the big noise
var_1 = 1e-6; rho_corrupt = .1;
% seeded decoding matrix parameters
alpha = 0.205; L = 10; J = .2; alpha1 = .22; alpha2 = alpha - alpha1 ./ L; W = 3;
% homogeneous decoding matrix parameters
alpha = 0.3;
% do you want to use l1magic pack (that must be added to the path)
% or the inluded l1 AMP solver?
% Default is zero, Be sure to download it and put if you use it: 
% http://users.ece.gatech.edu/~justin/l1magic/
l1_magic = 0; l1_AMP = 1 - l1_magic;

%% Initialisation
M = ceil(N ./ (1 - alpha) ); phant = phantom(sqrt(N) );
I_L1 = zeros(size(phant) ); I_CS = zeros(size(phant) ); I_id = zeros(size(phant) );
% Generate the measurement matrix and the coding one
if (seeded == 1)
    [F, Jmat, Mvec, Nvec, alpha] = Create_MAT1D_gauss_Wplus(M, L, J, alpha1, alpha2, W); F = full(F)';
else
    F = randn(M,M - N);
end
A = null(F');

disp('phantom decoding')
disp('gamma')
disp(1./ (1 - alpha))

%% Generate the signal and the perfect measure
S = phant;
rp = randperm(N); S = reshape(S, N, 1); S_mixed = S(rp);
comp = 0; while (min(size(A) ) > max(size(S_mixed) ) ); S_mixed = [S_mixed; randn]; comp = comp + 1; end;
Y_perfect = A * S_mixed; MM = max(size(Y_perfect) );

%% Error due to noisy channel transmission
k = floor(rho_corrupt .* MM);
rp2 = randperm(MM); error = [[randn(k, 1); zeros(MM - k, 1)] ]; error_ = error(rp2); error = error_;
small_error = randn(MM, 1) .* sqrt(var_1);
Y_err = Y_perfect + small_error + error;
Y_null = F' * Y_err;

%% Ideal coding and decoding
Y_id = Y_perfect + small_error;
X_id = A \ Y_id; X_id_(rp) = X_id(1:max(size(X_id) ) - comp, 1); X_id = X_id_;

%% Naive decoding (with no error correction)
X_naive = A \ Y_err; X_naive_(rp) = X_naive(1:max(size(X_naive) ) - comp, 1); X_naive = X_naive_;

%% BPCS coding and decoding by reconstruction of the sparse error vector
% algorithm properties
My = CSBP_Solver_Opt();
My.nb_iter = 1000;
My.print = 50;
My.conv = 10^(-8);
My.signal_rho = rho_corrupt;
My.signal = error + small_error;
My.dump_mes = 0.5;
My.method = 'AMP';
% 2Gauss
My.m_2_gauss = 0;
My.var_1_gauss = var_1;
My.var_2_gauss = 1;

%% AMP decoding with 2Gauss prior
My.prior = '2Gauss';
[results2G, n_and_e] = CSBP_Solver(Y_null', F', My);
Y_CS = Y_err - results2G.av_mess';
X_CS = A \ Y_CS; X_CS_(rp) = X_CS(1:max(size(X_CS) ) - comp, 1);  X_CS = X_CS_;

%% L1
if (l1_AMP == 1)
    My.prior = 'L1';
    [resultsL1 ,n_and_e] = CSBP_Solver(Y_null', F', My);
    Y_L1 = Y_err - resultsL1.av_mess';
    X_L1 = A \ Y_L1; X_L1_(rp) = X_L1(1 : max(size(X_L1) ) - comp, 1); X_L1 = X_L1_;
else
    % Dantzig selector
    err_L1 = l1dantzig_pd(F * Y_null, F', [], Y_null, 3e-3, 1e-3,1000);
    % regression step
    err_index = abs(err_L1) > sqrt(var_1 .* M); F_t = F';
    F_index_t = F_t(:,err_index);
    err_L1 = zeros(size(Y_err) ); err_L1(err_index) = F_index_t \ Y_null;
    Y_L1 = Y_err - err_L1;
    X_L1 = A \ Y_L1; X_L1_(rp) = X_L1(1 : max(size(X_L1) ) - comp, 1); X_L1 = X_L1_;
end

%% Results
S_unmixed(rp) = S(1 : max(size(S_mixed) ) - comp, 1); S = S_unmixed;
rho_id_CS = norm(X_CS - S, 2) ./ norm(X_id - S, 2);
rho_id_L1 = norm(X_L1 - S, 2) ./ norm(X_id - S, 2);
mse_id = mean((X_id - S).^2);
mse_CS = mean((X_CS - S).^2);
mse_L1 = mean((X_L1 - S).^2);

results_(1,1) = {reshape(X_id, sqrt(N), sqrt(N) )};
results_(2,1) = {reshape(X_CS, sqrt(N), sqrt(N) )};
results_(3,1) = {reshape(X_L1, sqrt(N), sqrt(N) )};
results_(4,1) = {[mse_id, mse_CS, mse_L1, rho_id_CS, rho_id_L1]};
results_(5,1) = {[max(size(S) ), MM, k, min(size(F) ) ./ MM, rho_corrupt, L, J, alpha1, alpha2, W]};

I_CS = reshape(X_CS, sqrt(N), sqrt(N) );
I_L1 = reshape(X_L1, sqrt(N), sqrt(N) );
I_id = reshape(X_id, sqrt(N), sqrt(N) );
I_naive = reshape(X_naive, sqrt(N), sqrt(N) );

%% Plots
clf
subplot(3, 2, 1); imshow(phant); xlabel('Original Image');
subplot(3, 2, 2); imshow(I_id); xlabel('Ideal');
subplot(3, 2, 3); imshow(I_CS); xlabel('AMP 2Gauss');
subplot(3, 2, 4); imshow(I_L1); xlabel('L1');
subplot(3, 2, 5); imshow(I_naive); xlabel('No error correction');

%% Save
save('phantom_decoding','I_id','I_CS','I_L1','I_naive','results_');