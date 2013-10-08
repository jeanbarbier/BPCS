alphaMin = 0.25;
alphaMax = 0.95;
dAlpha = 0.05;
rho = 0.1;
varNoise = 0;
convDE = 1e-10;
rep = 0; cav = 1;
step = 0.1; % the step to increase the lambda in the runs of the algo
average = 0;
variance = 1;
N = 1e4;

% algorithm properties
My = CSBP_Solver_Opt();
My.nb_iter = 100;
My.save_speed = 1;
My.save_memory = 1 - My.save_speed;
My.print = 10;
My.conv = 1e-15;
My.signal_rho = rho;
My.var_noise = 0;
My.dump_mes = 0.5;
My.prior = 'SparseGauss';
My.option_noise = 0;
My.remove_mean = 0;
My.method = 'AMP';
My.m_gauss = average;
My.var_gauss = variance;

t = 1; iMin = -5;
for (alpha = alphaMin : dAlpha : alphaMax)
    
    E = 0.1; V = E; DE = 1 + convDE; lambda = 0;
    while (DE > convDE);
        
        % Density evolution
        [E_new_rep,V_new_rep,E_new_cav,V_new_cav] = DensEvoSparseGauss(rep,cav,rho,alpha,varNoise,E,V);
        
        % Compute the maximum lambda (mean of the measurement matrix) above which there is an instability
        lambda = max(lambda, abs((varNoise + E) ./ (alpha .* E_new_cav) ) );
        
        DE = abs(E_new_cav - E); E = E_new_cav;
        
    end
    disp('lambda instability density evolution treshold');
    disp(lambda);
    lambdaDE(t) = lambda;
    
    conv = 1; i = iMin; ii = 0;
    while (conv == 1);
        
        meanG = (lambda + i .* step) ./ N;
        
        S = S_SparseGauss(N, rho, average, variance);
        G = randn(ceil(alpha .* N), N) ./ sqrt(N) + meanG;
        Y = G * S';
        
        My.signal = S;
        
        % reconstruction
        [results, n_and_e] = CSBP_Solver(Y, G, My);
        
        if (n_and_e.true_error > 1e-5);
            if (ii > 0); conv = 0; lambdaAlgo = lambda + i .* step;
            else i = i + iMin - 1;
            end
        end
        i = i + 1;
        ii = ii + 1;
        
    end
    
    disp('lambda algorithm instability treshold');
    disp(lambdaAlgo);
    lambdaAL(t) = lambdaAlgo;
    
    alpha_(t) = alpha;
    t = t + 1;
    
end