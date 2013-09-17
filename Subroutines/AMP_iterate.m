if (opt.save_speed == 1)
    if (opt.remove_mean == 0)
        V_new = (G2 * prior.var_mess')';
        W_new = dumping(W,(G * prior.av_mess')' - ((Y - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        var_1 = (1 ./ (n_and_e.var_noise + V_new) ) * G2;
        var_2 = ((Y - W_new) ./ (n_and_e.var_noise + V_new) ) * G;
    elseif (opt.remove_mean == 1)
        V_new = (G2 * prior.var_mess' + Gm1d2 * prior.var_mess' - 2 .* GGm * prior.var_mess')';
        W_new = dumping(W,(G * prior.av_mess' - Gmean1d * prior.av_mess')' - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
        PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
        var_1 = max(1e-18,PRELI * G2 + sum(PRELI) * Gm1d2 - 2 .* PRELI * GGm);
        var_2 = PRELI2 * G - sum(PRELI2) * Gmean1d;
    else
        V_new = (G2 * prior.var_mess' + Gm2 .* sum(prior.var_mess) - (2 .* GGm) * prior.var_mess')';
        W_new = dumping(W,(G * prior.av_mess' - Gmean .* sum(prior.av_mess))' - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
        PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
        var_1 = PRELI * G2 + sum(PRELI) .* Gm2 - 2 .* PRELI * GGm;
        var_2 = PRELI2 * G - sum(PRELI2) .* Gmean;
    end
end

if (opt.save_memory == 1)
    if (opt.remove_mean == 0)
        V_new = (G.^2 * prior.var_mess')';
        W_new = dumping(W,(G * prior.av_mess')' - ((Y - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        var_1 = (1 ./ (n_and_e.var_noise + V_new) ) * G.^2;
        var_2 = ((Y - W_new) ./ (n_and_e.var_noise + V_new) ) * G;
    elseif (opt.remove_mean == 1)
        V_new = (G.^2 * prior.var_mess' + Gm1d2 * prior.var_mess' - (2 .* G .* (ones(M,1) * sum(G) ./ M) ) * prior.var_mess')';
        W_new = dumping(W,(G * prior.av_mess' - Gmean1d * prior.av_mess')' - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
        PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
        var_1 = max(1e-18,PRELI * G.^2 + sum(PRELI) * Gm1d2 - 2 .* PRELI * (G .* (ones(M,1) * sum(G) ./ M) ) );
        var_2 = PRELI2 * G - sum(PRELI2) * Gmean1d;
    else
        V_new = (G.^2 * prior.var_mess' + Gm2 .* sum(prior.var_mess) - (2 .* G .* Gmean) * prior.var_mess')';
        W_new = dumping(W,(G * prior.av_mess' - Gmean .* sum(prior.av_mess))' - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
        PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
        var_1 = PRELI * G.^2 + sum(PRELI) * Gm2 - 2 .* PRELI * (G .* Gmean);
        var_2 = PRELI2 * G - sum(PRELI2) .* Gmean;
    end
end

if (opt.learn_d == 1); prior_temp.R = var_2 ./ var_1 + prior_temp.av_mess; prior_temp.S2 = 1 ./ var_1; prior_temp = F(prior_temp);
else prior.R = var_2 ./ var_1 + prior.av_mess; prior.S2 = 1 ./ var_1; prior = F(prior); V = V_new; W = W_new; end;
% Iteration done ----