if (opt.print > 0)
    switch opt.prior
        case 'SparseGauss'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m  var   convergence')
            else
                disp('iter   rho   noise   m  var   convergence   error')
            end
        case 'SparseGaussPositive'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m  var   convergence')
            else
                disp('iter   rho   noise   m  var   convergence   error')
            end
        case '2Gauss'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m_1   m_2   var_1   var_2   convergence')
            else
                disp('iter   rho   noise   m_1   m_2   var_1   var_2   convergence   error')
            end
        case 'SparseExponential'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   expo   convergence')
            else
                disp('iter   rho   noise   expo   convergence   error')
            end
        case 'SparseConstant'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   c_down   c_up   convergence')
            else
                disp('iter   rho   noise   c_down   c_up   convergence   error')
            end
        case 'SparseBinary'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   convergence')
            else
                disp('iter   rho   noise   convergence   error')
            end
        case 'L1'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   convergence')
            else
                disp('iter   rho   noise   convergence   error')
            end
        case 'Laplace'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   beta   convergence')
            else
                disp('iter   rho   noise   beta   convergence   error')
            end
        case 'Binary1'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   convergence')
            else
                disp('iter   rho   noise   convergence   error')
            end
        otherwise
            disp('unknown prior')
    end
end