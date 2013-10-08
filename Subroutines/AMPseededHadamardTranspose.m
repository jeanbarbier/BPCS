Gvar = MultSeededHadamardTranspose2(prior.var_mess, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError)';
Gav = MultSeededHadamardTranspose2(prior.av_mess, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError)';

if (opt.remove_mean == 0)
    V_new = MultSeededHadamardTransposeSquarred2(prior.var_mess, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock)';
    W_new = dumping(W,Gav - ((Y - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
    V_new = dumping(V,V_new,opt.dump_mes);
    var_1 = MultSeededHadamardSquarred2((1 ./ (n_and_e.var_noise + V_new) ), opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock)';
    var_2 = MultSeededHadamard2(((Y - W_new) ./ (n_and_e.var_noise + V_new) ), opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError)';
elseif (opt.remove_mean > 0)
    GvarGmean = zeros(1, opt.M);
    partVnew = GvarGmean;
    partWnew = GvarGmean;
    xx = 1;
    for l = 1 : opt.numBlockL
        a = MultSeededHadamardTranspose2(prior.var_mess .* Gmean(l, :), opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError)';
        GvarGmean(xx : xx + opt.Mblock(l) - 1) = a(xx : xx + opt.Mblock(l) - 1);
        partVnew(xx : xx + opt.Mblock(l) - 1) = Gmean2(l, :) * prior.var_mess' - 2 * GvarGmean(xx : xx + opt.Mblock(l) - 1);
        partWnew(xx : xx + opt.Mblock(l) - 1) = Gmean(l, :) * prior.av_mess';
        xx = xx + opt.Mblock(l);
    end
    
    V_new = MultSeededHadamardTransposeSquarred2(prior.var_mess, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock)' + partVnew;
    W_new = dumping(W,(Gav - partWnew) - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
    V_new = dumping(V,V_new,opt.dump_mes);
    
    PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
    PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
    
    part1var_1 = zeros(1, opt.N);
    part2var_1 = part1var_1;
    part1var_2 = part1var_1;
    xx = 1;
    for l = 1 : opt.numBlockL
        PRELIbyBlock = sum(PRELI(xx : xx + opt.Mblock(l) - 1) );
        PRELI2byBlock = sum(PRELI2(xx : xx + opt.Mblock(l) - 1) );
        part1var_1 = part1var_1 + Gmean2(l, :) .* PRELIbyBlock;
        PRELIzero = zeros(1, opt.M);
        PRELIzero(xx : xx + opt.Mblock(l) - 1) = PRELI(xx : xx + opt.Mblock(l) - 1);
        u = MultSeededHadamard2(PRELIzero, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError)';
        part2var_1 = part2var_1 + Gmean(l, :) .* u;
        part1var_2 = part1var_2 + Gmean(l, :) .* PRELI2byBlock;
        xx = xx + opt.Mblock(l);
    end
    
    var_1 = max(1e-18, part1var_1 + MultSeededHadamardSquarred2(PRELI, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock)' - 2 * part2var_1);
    var_2 = MultSeededHadamard2(PRELI2, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError)' - part1var_2;
end

prior.R = var_2 ./ var_1 + prior.av_mess; prior.S2 = 1 ./ var_1; prior = F(prior); V = V_new; W = W_new;
% Iteration done ----