% Printing to screen
switch opt.prior
    
    case ('SparseGauss')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('SparseGaussPositive')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('SparseConstant')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('2Gauss')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 prior.param_3 prior.param_4 n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 prior.param_3 prior.param_4 n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('SparseExponential')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('Laplace')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('SparseBinary')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('L1')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('Binary1')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise n_and_e.convergence n_and_e.true_error] ) );
        end
    case ('Complex')
        if (max(size(opt.signal) ) < 2)
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence] ) );
        else
            PR = sprintf('%d %e %e %e %e %e %e %e',full([t prior.rho n_and_e.var_noise prior.param_1 prior.param_2 n_and_e.convergence n_and_e.true_error] ) );
        end
    otherwise
        disp('unknown prior')
end

disp(PR)