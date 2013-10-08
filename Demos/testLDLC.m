clear all

N = 100;
connec = 3;
step = 0.3;


% capa = exp(-log(2 .* pi) + log(abs(det(G) ) ) .* 2 ./ N - 1);


compute = 1; capa = 1 ./ (2 .* pi .* exp(1) ); varNoise = capa; i = 1;
while (compute == 1)
    
    G = zeros(N);
    for i = 1 : N;
        rp = randperm(N);
        G(i, rp(1 : connec) ) = sign(randn(1, connec) );
    end
    G = G ./ det(G);
    S = S_SparseBinary(N, 0.5);
    S = S - (S == 0);
    
    Y = G * S' + sqrt(varNoise) .* randn(min(size(G) ), 1);
    
    % algorithm properties
    My = CSBP_Solver_Opt();
    My.nb_iter = 1000;
    My.save_speed = 1;
    My.save_memory = 1 - My.save_speed;
    My.print = 200;
    My.learn = 0;
    My.signal_rho = 1;
    My.var_noise = varNoise;
    My.signal = S;
    My.dump_learn = 0.5;
    My.dump_mes = 0.5 ;
    My.prior = 'Binary1';
    My.option_noise = 0;
    My.remove_mean = 0;
    My.method = 'AMP';
    My.conv = 0;
    
    % reconstruction
    [results,n_and_e] = CSBP_Solver(Y,G,My);
    
    varNoise = varNoise .* exp(-step);
    if (abs(log(varNoise ./ capa) ./ log(10) ) > 3); compute = 0; end
    error(i) = 1 ./ N .* sum(abs(results.av_mess - S) > 1e-12);
    var(i) = varNoise;
    dist(i) = log(capa ./ varNoise) ./ log(10);
    i = i + 1;
    
    log(capa ./ varNoise) ./ log(10)
    %     capa ./ varNoise
    1 ./ N .* sum(abs(results.av_mess - S) > 0)
end