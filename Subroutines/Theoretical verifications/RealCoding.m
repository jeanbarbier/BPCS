clear all
BG_ = 1; n_iter = 1e5; N = 16; rho_corrupt = .1; tt = 1;

%% 2Gauss
if (BG_ == 1)
    
    for gamma = 1.6:-0.025:1.2
        
        alpha = 1 - 1 ./ gamma; M = ceil(N ./ (1 - alpha) ); k = floor(rho_corrupt .* M); true_gamma = M / N;
        t = 1;
        while (t <= n_iter)
            
            %% Coding matrix
            G_null = randn(M,M - N); G = null(G_null'); S = randn(N,1);
            
            %% error
            error = intrlv([randn(k,1); zeros(M - k,1)], randperm(M) );
            
            %% Generate the signal, the perfect measure and the small error
            S = randn(N,1);
            Y_perfect = G * S;
            %             var_1 = (median(Y_perfect) ./ 16).^2;
            var_1 = 1e-6;
            small_error = randn(M,1) .* sqrt(var_1);
            
            %% Ideal coding and decoding
            %             Y_id = Y_perfect + small_error;
            %             X_id = G \ Y_id;
            
            %% Oracle coding and decoding
            %             j = 1; G_or = [];
            %             for i = 1 : M
            %                 if (error(i,1) == 0); Y_or(j,1) = Y_id(i,1); G_or = [G_or; G(i,:)]; j = j + 1; end;
            %             end
            %             X_or = G_or \ Y_or;
            
            %% BPCS coding and decoding by reconstruction of the sparse error vector
            Y_err = Y_perfect + small_error + error;
            Y_null = G_null' * Y_err;
            
            % algorithm properties
            My = CSBP_Solver_Opt();
            My.nb_iter = 1000;
            My.save_speed = 1;
            My.save_memory = 1 - My.save_speed;
            My.print = -1;
            My.conv = 10^(-10);
            My.learn = 0;
            My.learn_d = 0;
            My.signal_rho = rho_corrupt; % (-1 -> measure_rate ./ 10)
            My.calta = 0;
            My.var_noise = 0;
            My.signal = error + small_error;
            My.dump_learn = 0.;
            My.dump_mes = 0.5;
            My.prior = '2Gauss';
            My.option_noise = 0;
            My.remove_mean = 0;
            My.method = 'AMP';
            % 2Gauss
            My.m_2_gauss = 0;
            My.var_1_gauss = var_1;
            My.var_2_gauss = 1 + var_1;
            
            [results,n_and_e] = CSBP_Solver(Y_null',G_null',My);
            Y_CS = Y_err - results.av_mess';
            X_CS = G \ Y_CS;
            
            %% Results
            %             rho_id = norm(X_CS - S,2) ./ norm(X_id - S,2);
            %             rho_or = norm(X_CS - S,2) ./ norm(X_or - S,2);
            %             mse_id = mean((X_id - S).^2);
            %             mse_or = mean((X_or - S).^2);
            mse_CS = mean((X_CS - S).^2);
            
            %% Display and save results
            %             MSE__id(1,t) = mse_id;
            %             MSE__or(1,t) = mse_or;
            MSE__CS(1,t) = mse_CS;
            %             rho__id(1,t) = rho_id;
            %             rho__or(1,t) = rho_or;
            %             var_small_error(1,t) = var_1;
            
            [gamma, t]
            
            [statcode, result] = system(['ls stop',num2str(N)]);
            if statcode == 0
                results_2G(1,tt) = {MSE__CS'}; results_2G(2,tt) = {[true_gamma, gamma, k, rho_corrupt, var_1]};
                keyboard
                system(['rm stop',num2str(N)]);
            end
            t = t + 1;
        end
        
        %         results_2G(1,tt) = {[MSE__id',MSE__or',MSE__CS',rho__id',rho__or']}; results_2G(2,tt) = {[N, M, k, alpha, rho_corrupt]}; results_2G(3,tt) = {var_small_error};
        %         results_2G(1,tt) = {[MSE__id',MSE__CS',rho__id']}; results_2G(2,tt) = {[N, M, k, alpha, rho_corrupt]}; results_2G(3,tt) = {var_small_error};
        results_2G(1,tt) = {MSE__CS'}; results_2G(2,tt) = {[true_gamma, gamma, k, rho_corrupt, var_1]};
        M = ceil(N ./ (1 - alpha) ); k = floor(rho_corrupt .* M); tt = tt + 1;
        
    end
    
    save(['Full2G',num2str(N)],'results2G')
    
end

% clear all
% L1_ = 0; n_iter = 500; N_max = 256; alpha = .5; rho_corrupt = .1;
%
% %% L1
% if (L1_ == 1)
%     tt = 1; N = 256; M = ceil(N ./ (1 - alpha) ); k = floor(rho_corrupt .* M);
%
%     while (N <= N_max)
%
%         for t = 1 : n_iter
%
%             %% Coding matrix
%             G_null = randn(M,M - N); G = null(G_null'); S = randn(N,1);
%
%             %% error
%             error = intrlv([randn(k,1); zeros(M - k,1)], randperm(M) );
%
%             %% Generate the signal, the perfect measure and the small error
%             Y_perfect = G * S;
%             %             var_1 = (median(Y_perfect) ./ 16).^2;
%             var_1 = 1e-6;
%             small_error = randn(M,1) .* sqrt(var_1);
%
%             %% Ideal coding and decoding
%             Y_id = Y_perfect + small_error;
%             X_id = G \ Y_id;
%
%             %% Oracle coding and decoding
%             j = 1; G_or = [];
%             for i = 1 : M
%                 if (error(i,1) == 0); Y_or(j,1) = Y_id(i,1); G_or = [G_or; G(i,:)]; j = j + 1; end;
%             end
%             X_or = G_or \ Y_or;
%
%             %% BPCS coding and decoding by reconstruction of the sparse error vector
%             Y_err = Y_perfect + error + small_error;
%             Y_null = G_null' * Y_err;
%
%             %             % algorithm properties
%             %             My = CSBP_Solver_Opt();
%             %             My.nb_iter = 1000;
%             %             My.save_speed = 1;
%             %             My.save_memory = 1 - My.save_speed;
%             %             My.print = 100;
%             %             My.conv = 10^(-20);
%             %             My.learn = 0;
%             %             My.learn_d = 0;
%             %             My.signal_rho = rho_corrupt; % (-1 -> measure_rate ./ 10)
%             %             My.calta = 0;
%             %             My.var_noise = 0;
%             %             My.signal = error + small_error;
%             %             My.dump_learn = 0.;
%             %             My.dump_mes = 0.5;
%             %             My.prior = 'SparseGauss';
%             %             My.option_noise = 0;
%             %             My.remove_mean = 0;
%             %             My.method = 'AMP';
%             %             % L1
%             %             My.min = -10;
%             %             My.max = 10;
%             %             % 2Gauss
%             %             My.m_2_gauss = 0;
%             %             My.var_1_gauss = var_1;
%             %             My.var_2_gauss = 1;
%             %
%             %             [results,n_and_e] = CSBP_Solver(Y_null',G_null',My);
%             %
%             %             Y_CS = Y_err - results.av_mess';
%             %             X_CS = G \ Y_CS;
%             %
%             %             % regression step
%             %             err_index = abs(results.av_mess') > sqrt(var_1); G_null_t = G_null';
%             %             G_index_t_null = G_null_t(:,err_index);
%             %             err_CS = zeros(size(Y_err) );
%             %             err_CS(err_index) = G_index_t_null \ Y_null;
%             %             Y_CS = Y_err - err_CS;
%             %             X_CS = G \ Y_CS;
%
%             %% L1 magic
%             %             % L1 optimization of the error
%             %             disp('L1_1')
%             %             err_L1_1 = l1eq_pd(G_null * Y_null, G_null', [], Y_null, 1e-3, 1000);
%             %             size(Y_null)
%             %             Y_L1_1 = Y_err - err_L1_1;
%             %             X_L1_1 = G \ Y_L1_1;
%             %
%             %         % Error optimization
%             %         disp('L1_2')
%             %         X_L1_2 = l1decode_pd(G' * Y_err, G, [], Y_err, 1e-4, 1000);
%             %
%
%             % Dantzig selector
%             Y = Y_null; F = G_null';
%             err_L1_3 = l1dantzig_pd(F' * Y, F, [], Y, 3e-3, 1e-3,1000);
%             % regression step
%             err_index = abs(err_L1_3) > sqrt(var_1 .* M);
%             F_index = F(:,err_index);
%             err_L1_3 = zeros(size(Y_err) );
%             err_L1_3(err_index) = F_index \ Y;
%             Y_L1_3 = Y_err - err_L1_3;
%             X_L1_3 = G \ Y_L1_3;
%             %
%             %         % Quad constraint L1 optimization of error
%             %         disp('L1_4')
%             %         err_L1_4 = l1qc_logbarrier(G_null * Y_null, G_null', [], Y_null, 3e-3);
%             %         Y_L1_4 = Y_err - err_L1_4;
%             %         X_L1_4 = G \ Y_L1_4;
%
%             %% Results
%             rho_id = norm(X_CS - S,2) ./ norm(X_id - S,2);
%             rho_or = norm(X_CS - S,2) ./ norm(X_or - S,2);
%             mse_id = mean((X_id - S).^2);
%             mse_or = mean((X_or - S).^2);
%             mse_CS = mean((X_CS - S).^2);
%
%             %% Display and save results
%             MSE__id(1,t) = mse_id;
%             MSE__or(1,t) = mse_or;
%             MSE__CS(1,t) = mse_CS;
%             rho__id(1,t) = rho_id;
%             rho__or(1,t) = rho_or;
%             var_small_error(1,t) = var_1;
%
%             mse_CS
%             t
%
%         end
%
%         results_L1(1,tt) = {[MSE__id',MSE__or',MSE__CS',rho__id',rho__or']}; results_L1(2,tt) = {[N, M, k, alpha, rho_corrupt]}; results_L1(3,tt) = {var_small_error};
%         %         results_L1(1,tt) = {[MSE__CS']}; results_L1(2,tt) = {[N, M, k, alpha, rho_corrupt]}; results_L1(3,tt) = {var_small_error};
%         N = N .* 2; M = ceil(N ./ (1 - alpha) ); k = floor(rho_corrupt .* M); tt = tt + 1;
%
%     end
%
%     save('FullL1','results_L1')
%
% end