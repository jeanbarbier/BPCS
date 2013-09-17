clear all
N = 64.^2; var_1 = 1e-6; L = 10; J = .2; alpha1 = .22; alpha2 = .205 - .022; alpha = alpha1 ./ L + alpha2; W = 3; rho_corrupt = .1; M = ceil(N ./ (1 - alpha) );
Lena = imread('Lena256.png');
I_L1 = zeros(size(Lena) ); I_CS = zeros(size(Lena) ); I_id = zeros(size(Lena) );

% Generate the measurement matrix and the coding one
[G_null, Jmat, Mvec, Nvec] = Create_MAT1D_gauss_Wplus(M, L, J, alpha1, alpha2, W);
G_null = full(G_null)'; G = null(G_null');

disp('ok')

for ll = 1 : max(size(Lena) ) ./ sqrt(N)
    for cc = 1 : max(size(Lena) ) ./ sqrt(N)
        
        %% Generate the  signal, the perfect measure and the small error
        S = double(Lena((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1).* sqrt(N) + 1: cc .* sqrt(N) ) ) ./ 256;
        rp = randperm(N); S = intrlv(reshape(S, N, 1), rp);
        comp = 0; while (min(size(G) ) > max(size(S) ) ); S = [S; randn]; comp = comp + 1; end;
        Y_perfect = G * S; MM = max(size(Y_perfect) );
        small_error = randn(MM, 1) .* sqrt(var_1);
        
        %% Ideal coding and decoding
        Y_id = Y_perfect + small_error;
        X_id = G \ Y_id; X_id = deintrlv(X_id(1:max(size(X_id) ) - comp, 1), rp);
        
        %% error
        k = floor(rho_corrupt .* MM);
        error = intrlv([randn(k, 1); zeros(MM - k, 1)], randperm(MM) );
        
        %% BPCS coding and decoding by reconstruction of the sparse error vector
        Y_err = Y_perfect + small_error + error;
        Y_null = G_null' * Y_err;
        
        % algorithm properties
        My = CSBP_Solver_Opt();
        My.nb_iter = 1000;
        My.save_speed = 1;
        My.save_memory = 1 - My.save_speed;
        My.print = 50;
        My.conv = 10^(2);
        My.learn = 0;
        My.learn_d = 0;
        My.signal_rho = rho_corrupt; % (-1 -> measure_rate ./ 10)
        My.calta = 0;
        My.var_noise = 0;
        My.signal = error;
        My.dump_learn = 0.;
        My.dump_mes = 0.5;
        My.option_noise = 0;
        My.remove_mean = 0;
        My.method = 'AMP';
        % 2Gauss
        My.m_2_gauss = 0;
        My.var_1_gauss = var_1;
        My.var_2_gauss = 1;
        % L1
        My.min = -10;
        My.max = 10;
        
        %% 2Gauss
        My.prior = '2Gauss'; err_CS = 1;
        while (err_CS >= 1e-4)
            [G_null, Jmat, Mvec, Nvec] = Create_MAT1D_gauss_Wplus(M, L, J, alpha1, alpha2, W);
            G_null = full(G_null)'; G = null(G_null'); Y_null = G_null' * Y_err;
            [results2G, n_and_e] = CSBP_Solver(Y_null', G_null', My); err_CS = n_and_e.true_error; fail = fail + 1;
        end
        Y_CS = Y_err - results2G.av_mess';
        X_CS = (G \ Y_CS); X_CS = deintrlv(X_CS(1 : max(size(X_CS) ) - comp, 1), rp);
        
        %% L1
        %         My.prior = 'L1';
        %         [resultsL1,n_and_e] = CSBP_Solver(Y_null', G_null', My);
        %         Y_L1 = Y_err - resultsL1.av_mess';
        %         X_L1 = (G \ Y_L1); X_L1 = deintrlv(X_L1(1 : max(size(X_L1) ) - comp, 1), rp);
        
        %         err_L1 = l1eq_pd(G_null * Y_null, G_null', [], Y_null, 1e-3, 1000);
        %         Y_L1 = Y_err - err_L1;
        %         X_L1 = G \ Y_L1; X_L1 = deintrlv(X_L1(1 : max(size(X_L1) ) - comp, 1), rp);
        
        % Dantzig selector
        err_L1_3 = l1dantzig_pd(G_null * Y_null, G_null', [], Y_null, 3e-3, 1e-3,1000);
        % regression step
        err_index = abs(err_L1_3) > sqrt(var_1 .* M); G_null_t = G_null';
        G_index_t_null = G_null_t(:,err_index);
        err_L1_3 = zeros(size(Y_err) );
        err_L1_3(err_index) = G_index_t_null \ Y_null;
        mean((err_L1_3 - (2 .* Y_perfect .* error) ).^2)
        Y_L1_3= Y_err - err_L1_3;
        X_L1 = G \ Y_L1_3; X_L1 = deintrlv(X_L1(1 : max(size(X_L1) ) - comp, 1), rp);
        
        %% Results
        S = deintrlv(S(1 : max(size(S) ) - comp, 1), rp);
        rho_id_CS = norm(X_CS - S, 2) ./ norm(X_id - S, 2);
        rho_id_L1 = norm(X_L1 - S, 2) ./ norm(X_id - S, 2);
        mse_id = mean((X_id - S).^2);
        mse_CS = mean((X_CS - S).^2);
        mse_L1 = mean((X_L1 - S).^2);
        
        results_(1, ll, cc) = {reshape(X_id, sqrt(N), sqrt(N) )};
        results_(2, ll, cc) = {reshape(X_CS, sqrt(N), sqrt(N) )};
        results_(3, ll, cc) = {reshape(X_L1, sqrt(N), sqrt(N) )};
        results_(4, ll, cc) = {[mse_id, mse_CS, mse_L1, rho_id_CS, rho_id_L1]};
        results_(5, ll, cc) = {[max(size(S) ), MM, k, min(size(G_null) ) ./ MM, rho_corrupt, L, J, alpha1, alpha2, W]};
        
        I_CS((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1: cc .* sqrt(N) ) = reshape(X_CS, sqrt(N), sqrt(N) );
        I_L1((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1: cc .* sqrt(N) ) = reshape(X_L1, sqrt(N), sqrt(N) );
        I_id((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1: cc .* sqrt(N) ) = reshape(X_id, sqrt(N), sqrt(N) );
        
    end
end

%% Plots
% clf
% subplot(2, 2, 1); imshow(Lena); xlabel('Original Image');
% subplot(2, 2, 2); imshow(I_id); xlabel('Ideal');
% subplot(2, 2, 3); imshow(I_CS); xlabel('AMP 2Gauss');
% subplot(2, 2, 4); imshow(I_L1); xlabel('L1');

save('LENAseed','I_id','I_CS','I_L1','results_');