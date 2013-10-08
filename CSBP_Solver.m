function [prior, n_and_e, MSEt] = CSBP_Solver(Y, G, opt)
check_size();
set_parameters();
display_information();
initialisation();

% Construction of the prior dependent class, noise_and_error and dump ones
prior = Prior(opt.signal_rho, N, measure_rate, opt.learn, opt.prior, opt.dump_learn, R_init, S2_init, av_mess_init, var_mess_init, opt.method, prior_param{1}, prior_param{2}, prior_param{3}, prior_param{4}); F = str2func(prior.func); prior_temp = prior;
n_and_e = noise_and_error(opt.conv, opt.var_noise, opt.dump_learn);

% Starting main code
t = 0; print_to_screen(); t = 1;
while (t <= opt.nb_iter)
    
    switch (opt.method)
        case ('AMP'); AMP();
        case ('AMPtap'); AMPtap();
        case ('AMPtapB'); AMPtapB(); % ??
        case ('AMPseededHadamard'); AMPseededHadamard();
        case ('AMPseededHadamardTranspose'); AMPseededHadamardTranspose();
        case ('AMPseededHadamardTransposeA'); AMPseededHadamardTransposeA();
        case ('AMPseededFourier'); AMPseededFourier();
        case ('AMPcomplex'); AMPcomplex();
    end
    
    % Test of the convergence
    n_and_e = n_and_e.compute_convergence(prior.av_mess_old,prior.av_mess);
    if ((n_and_e.convergence < opt.conv) ); fprintf('Converged, convergence = %e',n_and_e.convergence); break; end;
    
    
    % Test of reconstruction on the fly knowing the original signal
    if (max(size(opt.signal) ) > 2)
        n_and_e = n_and_e.compute_true_MSE(opt.signal,prior.av_mess);
        if ((n_and_e.true_error < opt.conv) ); fprintf('Solution found, true error = %e',n_and_e.true_error); break; end;
    end
    
    % Learning of the noise if activated
    if (opt.option_noise == 1)
        n_and_e = n_and_e.learn_noise(Y,W_new,V_new);
        n_and_e.var_noise = dumping(n_and_e.var_noise_old, n_and_e.var_noise, opt.dump_learn);
        n_and_e.var_noise_old = n_and_e.var_noise;
    end
    
    if ((opt.print > 0) && (mod(t, opt.print) == 0) ); print_to_screen(); end
    
    % MSE by block
    if ((opt.MSEbyBlock > 0) && (mod(t, opt.MSEbyBlock) == 0) && (strcmp(opt.method, 'AMPseededHadamard') || strcmp(opt.method, 'AMPseededFourier') ) )
        MSE = MSEbyBloc(prior.av_mess, opt.signal, opt.numBlockC, opt.Nblock, opt.method);
        if(strcmp(opt.method, 'AMPhadamardSeeded') ); semilogy(MSE);
        else semilogy([1 : opt.numBlockC], MSE(1, :), 'r', [1 : opt.numBlockC], MSE(2, :), 'b'); end
        if (t == opt.MSEbyBlock); MSEmax = max(max(MSE) ); end
        axis([1, opt.numBlockC, 0, MSEmax] ); drawnow;
    end
    
    % MSE as a function of the iterations, to be compared with density evolution
    MSEt(t) = n_and_e.true_error;
    
    t = t + 1;
    
end

end