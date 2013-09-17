clear all
n_iter = 500; N_max = 4096; L = 10; J = .2; alpha1 = .22; alpha2 = 0.1630 - 0.022; alpha = alpha1 ./ L + alpha2; W = 3; rho_corrupt = .1; var_1 = 1e-6;

%% 2Gauss
tt = 1; N = 4096; M = ceil(N ./ (1 - alpha) );

%% Generate the coding matrix, the signal, the perfect measure and the small error
[G_null,Jmat,Mvec,Nvec,alpha_true] = Create_MAT1D_gauss_Wplus(M,L,J,alpha1,alpha2,W); G_null = full(G_null)'; G = null(G_null');

while (N <= N_max + 1)
    
    for t = 1 : n_iter
        
        %% Generate the signal
        S = randn(min(size(G) ),1);
        Y_perfect = G * S; MM = max(size(Y_perfect) );
        %         var_1 = (median(Y_perfect) ./ 16).^2;
        small_error = randn(MM,1) .* sqrt(var_1);
        
        %% Ideal coding and decoding
        %         Y_id = Y_perfect + small_error;
        %         X_id = G \ Y_id;
        
        %% error
        k = floor(rho_corrupt .* MM);
        error = intrlv([randn(k,1); zeros(MM - k,1)], randperm(MM) );
        
        %         %% Oracle coding and decoding
        %         j = 1; G_or = [];
        %         for i = 1 : MM
        %             if (error(i,1) == 0); Y_or(j,1) = Y_id(i,1); G_or = [G_or; G(i,:)]; j = j + 1; end;
        %         end
        %         X_or = G_or \ Y_or;
        
        %% BPCS coding and decoding by reconstruction of the sparse error vector
        Y_err = Y_perfect + small_error + error;
        Y_null = G_null' * Y_err;
        
        % algorithm properties
        My = CSBP_Solver_Opt();
        My.nb_iter = 1000;
        My.save_speed = 1;
        My.save_memory = 1 - My.save_speed;
        My.print = 100;
        My.conv = 10^(-10);
        My.learn = 0;
        My.learn_d = 0;
        My.signal_rho = rho_corrupt; % (-1 -> measure_rate ./ 10)
        My.calta = 0;
        My.var_noise = 0;
        My.signal = error;
        My.dump_learn = 0.;
        My.dump_mes = 0.5;
        My.prior = '2Gauss';
        My.option_noise = 0;
        My.remove_mean = 0;
        My.method = 'AMP';
        % 2Gauss
        My.m_2_gauss = 0;
        My.var_1_gauss = var_1;
        My.var_2_gauss = 1;
        
        [results,n_and_e] = CSBP_Solver(Y_null',G_null',My);
        Y_CS = Y_err - results.av_mess';
        X_CS = G \ Y_CS;
        
        %% Results
        %         rho_id = norm(X_CS - S,2) ./ norm(X_id - S,2);
        %         rho_or = norm(X_CS - S,2) ./ norm(X_or - S,2);
        %         mse_id = mean((X_id - S).^2);
        %         mse_or = mean((X_or - S).^2);
        mse_CS = mean((X_CS - S).^2);
        
        %% Display and save results
        %         rho__id(1,t) = rho_id;
        %         rho__or(1,t) = rho_or;
        %         MSE__id(1,t) = mse_id;
        %         MSE__or(1,t) = mse_or;
        MSE__CS(1,t) = mse_CS;
        
        var_small_error(1,t) = var_1;
        [alpha_true, t, mse_CS]
        
    end
    
    %     results_2G(1,tt) = {[MSE__id',MSE__or',MSE__CS',rho__id',rho__or']}; results_2G(2,tt) = {[max(size(S) ), MM, k, min(size(G_null) ) ./ MM, rho_corrupt, L, J, alpha1, alpha2, W]}; results_2G(3,tt) = {var_small_error};
    %     results_2G(1,tt) = {[MSE__id',MSE__CS',rho__id']}; results_2G(2,tt) = {[max(size(S) ), MM, k, min(size(G_null) ) ./ MM, rho_corrupt, L, J, alpha1, alpha2, W]}; results_2G(3,tt) = {var_small_error};
    results_2G(1,tt) = {[MSE__CS']}; results_2G(2,tt) = {[max(size(S) ), MM, k, min(size(G_null) ) ./ MM, rho_corrupt, L, J, alpha1, alpha2, W]}; results_2G(3,tt) = {var_small_error};
    N = N .* 2; M = ceil(N ./ (1 - alpha) ); tt = tt + 1;
    
end

save('2Gseed','results_2G')