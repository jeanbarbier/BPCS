clear all
alpha = 1e-2; dalpha = 5e-2; alpha_max = 1; t = 1;

while (alpha <= alpha_max)  
    % See option in generic_test_AlgoEvo
    generic_test_AlgoEvo();
    E_algo(1,t) = n_and_e.true_error;
    V_algo(1,t) = sum(results.var_mess) ./ size_S;
    alpha_(1,t) = alpha;
    t = t + 1; alpha = alpha + dalpha;
    
    if (n_and_e.true_error > 4);  t = t - 1; alpha = alpha - dalpha; end;
    
end
    