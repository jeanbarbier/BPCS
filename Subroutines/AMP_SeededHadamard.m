Gvar = MultSeededHadamard(prior.var_mess, opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock, opt.rpUp, opt.rpMiddle, opt.rpDown)';
Gav = MultSeededHadamard(prior.av_mess, opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock, opt.rpUp, opt.rpMiddle, opt.rpDown)';

if (opt.remove_mean == 0)
    V_new = Gvar;
    W_new = dumping(W,Gav - ((Y - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
    V_new = dumping(V,V_new,opt.dump_mes);
    var_1 = MultSeededHadamardTranspose((1 ./ (n_and_e.var_noise + V_new) ), opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock, opt.rpUp, opt.rpMiddle, opt.rpDown);
    var_2 = MultSeededHadamardTranspose(((Y - W_new) ./ (n_and_e.var_noise + V_new) ), opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock, opt.rpUp, opt.rpMiddle, opt.rpDown);
elseif (opt.remove_mean == 1)
    GvarGmean = MultSeededHadamard(prior.var_mess .* Gmean, opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock, opt.rpUp, opt.rpMiddle, opt.rpDown)';
    V_new = (MultSeededHadamardSquarred(prior.var_mess, opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock)' + Gm1d2 * prior.var_mess' - 2 .* GvarGmean);
    W_new = dumping(W,(Gav - Gmean1d * prior.av_mess') - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
    V_new = dumping(V,V_new,opt.dump_mes);
    PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
    PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
    var_1 = max(1e-18,sum(PRELI) * Gm1d2 + MultSeededHadamardTransposeSquarred(PRELI, opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock)' - 2 * Gmean .* MultSeededHadamardTranspose(PRELI, opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock, opt.rpUp, opt.rpMiddle, opt.rpDown)');
    var_2 = MultSeededHadamardTranspose(PRELI2, opt.Jup, opt.Jdown, opt.numBlock, opt.Mblock, opt.Nblock, opt.rpUp, opt.rpMiddle, opt.rpDown)' - sum(PRELI2) * Gmean1d;
else
    V_new = (Gvar' .* (1 - 2 .* Gmean) + Gm2 .* sum(prior.var_mess) );
    W_new = dumping(W,(Gvar' - Gmean .* sum(prior.av_mess)) - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
    V_new = dumping(V,V_new,opt.dump_mes);
    PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
    PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
    var_1 = Hadamard_multiply_vec_prime(PRELI,M,N,Index_random_part1) .* (1 - 2 .* Gmean) + sum(PRELI) .* Gm2;
    var_2 = Hadamard_multiply_vec_prime(PRELI2,M,N,Index_random_part1)  - sum(PRELI2) .* Gmean;
end

if (opt.learn_d == 1); prior_temp.R = var_2 ./ var_1 + prior_temp.av_mess; prior_temp.S2 = 1 ./ var_1; prior_temp = F(prior_temp);
else prior.R = var_2 ./ var_1 + prior.av_mess; prior.S2 = 1 ./ var_1; prior = F(prior); V = V_new; W = W_new; end;
% Iteration done ----