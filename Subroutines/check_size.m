% Simple checks on the sizes of Y and G and varG for bloc matrices
if ((strcmp(opt.method,'AMPseededHadamard') == 0) && (strcmp(opt.method,'AMPseededHadamardTranspose') == 0) && (strcmp(opt.method,'AMPseededHadamardTransposeA') == 0) && (strcmp(opt.method,'AMPseededFourier') == 0) )
    [M,N] = size(G);
    if ((M > N) && (alphaBig ~= 1) )
        G = G.';
        [M,N] = size(G);
    end
    
else
    M = opt.M; N = opt.N;
end

[a,b] = size(Y);
if (a > b)
    Y = Y.';
end
measure_rate = M ./ N;

if (strcmp(opt.method,'AMPhb') )
    [MM,NN] = size(opt.varG);
    if ((MM ~= M) || (NN ~= N) )
        if (((sum(opt.Nvec_bm) == []) > 0) | ((sum(opt.Mvec_bm) == []) > 0) ); error('varG has not the same size as G and thus opt.Mvec_bm and opt.Nvec_bm must be given as input');
        else opt.varG = convert2fullMN(opt.varG,opt.Mvec_bm,opt.Nvec_bm); end;
    end
    if (MM > NN)
        opt.varG = opt.varG.';
        [MM,NN] = size(opt.varG);
    end
end