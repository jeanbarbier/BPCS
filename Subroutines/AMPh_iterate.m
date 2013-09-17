if (opt.save_speed == 1)
    if (opt.remove_mean == 0)
        V_new = varG .* sum(prior.var_mess);
        W_new = dumping(W,(G * prior.av_mess')' - ((Y - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        S2_ = 1 ./ (sum(1 ./ (n_and_e.var_noise + V_new) ) .* M .* varG);
        R_p1 = (Y - W_new) ./ (n_and_e.var_noise + V_new) * G .* S2_;
    elseif (opt.remove_mean == 1)
        V_new = ((G2m + Gm1d2) * prior.var_mess' - 2 .* GGm * prior.var_mess')';
        W_new = dumping(W,(Geff * prior.av_mess')'- ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        S2_ = 1 ./ ((G2m + Gm1d2) .* sum(1 ./ (n_and_e.var_noise + V_new) ) - (2 ./ (n_and_e.var_noise + V_new) * GGm) );
        R_p1 = (((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ) * Geff) .* S2_;
    else
        V_new = ((G2m + Gm2) .* sum(prior.var_mess) - 2 .* GGm * prior.var_mess')';
        W_new = dumping(W,(Geff * prior.av_mess')' - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        S2_ = 1 ./ ((G2m + Gm2) .* sum(1 ./ (n_and_e.var_noise + V_new) ) - 2 ./ (n_and_e.var_noise + V_new) * GGm);
        R_p1 = (((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ) * Geff) .* S2_;
    end
end

if (opt.save_memory == 1)
    if (opt.remove_mean == 0)
        V_new = varG .* sum(prior.var_mess);
        W_new = dumping(W,(G * prior.av_mess')' - ((Y - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        S2_ = 1 ./ (sum(1 ./ (n_and_e.var_noise + V_new) ) .* M .* varG);
        R_p1 = (Y - W_new) ./ (n_and_e.var_noise + V_new)  * G .* S2_ ;
    elseif (opt.remove_mean == 1)
        V_new = ((G2m + Gm1d2) * prior.var_mess' - (2 .* G .* (ones(M,1) * sum(G) ./ M) ) * prior.var_mess')';
        W_new = dumping(W,(G * prior.av_mess' - Gmean1d * prior.av_mess')'- ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        S2_ = 1 ./ ((G2m + Gm1d2) .* sum(1 ./ (n_and_e.var_noise + V_new) ) - (2 ./ (n_and_e.var_noise + V_new) * (G .* (ones(M,1) * sum(G) ./ M) ) ) );
        R_p1 = (((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ) * G - Gmean1d .* sum((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ) ) .* S2_;
    else
        V_new = ((G2m + Gm2) .* sum(prior.var_mess) - (2 .* G .* (ones(M,1) * sum(G) ./ M) ) * prior.var_mess')';
        W_new = dumping(W,((G - Gmean) * prior.av_mess')' - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
        V_new = dumping(V,V_new,opt.dump_mes);
        S2_ = 1 ./ ((G2m + Gm2) .* sum(1 ./ (n_and_e.var_noise + V_new) ) - 2 ./ (n_and_e.var_noise + V_new) * (G .* (ones(M,1) * sum(G) ./ M) ) );
        R_p1 = (((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ) * (G - Gmean)) .* S2_;
    end
end

if (opt.learn_d == 1); prior_temp.R = R_p1 + prior_temp.av_mess; prior_temp.S2 = S2_; prior_temp = F(prior_temp);
else prior.R = R_p1 + prior.av_mess; prior.S2 = S2_; prior = F(prior); V = V_new; W = W_new; end;
% Iteration done ----