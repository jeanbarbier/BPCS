% Reading parameters
if (nargin <= 2); opt = CSBP_Solver_Opt(); end;

% Set parameters
Gmean = opt.Gmean; Ymean = opt.Ymean;
Nvec = opt.Nvec; if (Nvec == 1); Nvec = N; end;
Mvec = opt.Mvec; if (Mvec == 1); Mvec = M; end;

% Prior dependent quantities
switch opt.prior
    case 'SparseGauss'
        prior_param{1} = opt.m_gauss; prior_param{2} = opt.var_gauss; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'SparseGaussCut'
        prior_param{1} = opt.mC_gauss; prior_param{2} = opt.varC_gauss; prior_param{3} = opt.Cut; prior_param{4} = NaN;
    case 'SparseGaussPositive'
        prior_param{1} = opt.m_gaussP; prior_param{2} = opt.var_gaussP; prior_param{3} = NaN; prior_param{4} = NaN;
    case '2Gauss'
        prior_param{1} = opt.m_1_gauss; prior_param{2} = opt.m_2_gauss; prior_param{3} = opt.var_1_gauss; prior_param{4} = opt.var_2_gauss;
    case 'SparseExponential'
        prior_param{1} = opt.expo; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'SparseConstant'
        prior_param{1} = opt.c_down; prior_param{2} = opt.c_up; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'SparseBinary'
        prior_param{1} = NaN; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'L1'
        prior_param{1} = opt.min; prior_param{2} = opt.max; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'Laplace'
        prior_param{1} = opt.beta; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'Binary1'
        prior_param{1} = NaN; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    otherwise
        disp('unknown prior')
end

if (opt.signal_rho < 0); opt.signal_rho = measure_rate ./ 10; end;
[a,b] = size(opt.signal); if (a > b); opt.signal = opt.signal'; end;

if (strcmp(opt.method,'AMP') || strcmp(opt.method,'AMPh') )
    if (opt.save_speed == 1) % requires approximately 3 to 4 times the memory space needed to store the measurement matrix G but speed up the algorithm.
        if (opt.remove_mean == 0); varG = sum(sum(G.^2) ) ./ (M .* N); end;
        if (opt.remove_mean > 0)
            % Remove mean in G and Y
            if (max(size(Gmean) ) < 1)
                % Not provided by user, compute it!
                disp('Computing Mean');
                Gmean = sum(G) ./ M; % [1,N]
                Gmean1d = Gmean; % [1,N]
                if (opt.remove_mean == 2); Gmean = sum(Gmean) ./ N; else Gmean = ones(M,1) * Gmean; end;
            else
                disp('Using provided Means');
            end
            Yeff = Y - sum(Y) ./ M;
            Gm2 = Gmean.^2;
            Gm1d2 = Gmean1d.^2;
            GGm = G .* Gmean;
            Geff = G - Gmean;
        end
        if (strcmp(opt.method,'AMP') ); G2 = G.^2;
        elseif (strcmp(opt.method,'AMPh') && (opt.remove_mean == 2) ); G2m = sum(sum(G.^2) ) ./ (M .* N);
        elseif (strcmp(opt.method,'AMPh') && (opt.remove_mean == 1) ); G2m = sum(G.^2) ./ M;
        end;
    end
    
    if (opt.save_memory == 1) % requires 3 to 4 times less memory but the algorithm is slower because pre-computation are not allowed.
        if (opt.remove_mean == 0); varG = sum(sum(G.^2) ) ./ (M .* N); end;
        if (opt.remove_mean > 0)
            % Remove mean in G and Y
            if (max(size(Gmean) ) < 1)
                % Not provided by user, compute it!
                disp('Computing Mean');
                Gmean = sum(G) ./ M; % [1,N]
                Gmean1d = Gmean; % [1,N]
                if (opt.remove_mean == 2); Gmean = sum(Gmean) ./ N; end;
            else
                disp('Using provided Means');
            end
            Yeff = Y - sum(Y) ./ M;
            Gm2 = Gmean.^2;
            Gm1d2 = Gmean1d.^2;
        end
        if (strcmp(opt.method,'AMPh') && (opt.remove_mean == 2) ); G2m = sum(sum(G.^2) ) ./ (M .* N); varG = abs(G2m - Gm2);
        elseif (strcmp(opt.method,'AMPh') && (opt.remove_mean == 1) ); G2m = sum(G.^2) ./ M; varG = abs(G2m - Gm2);
        end;
    end
end


if ((opt.remove_mean == 0) && (opt.learn_d == 1) && (strcmp(opt.method,'AMP') || strcmp(opt.method,'AMPh') ) ); G2m = sum(G.^2) ./ M; end;

if (strcmp(opt.method,'AMPhadamard') == 1)
    load Index_random_part1;
    M = opt.M; N = opt.N;
    if (opt.remove_mean > 0); Gmean = 1 ./ M .* Hadamard_multiply_vec_prime(ones(1,M),M,N,Index_random_part1); end
    Gmean1d = Gmean; % [1,N]
    if (opt.remove_mean == 2); Gmean = sum(Gmean) ./ N; end
    Yeff = Y - sum(Y) ./ M;
    Gm2 = Gmean.^2;
    Gm1d2 = Gmean1d.^2;
end

if ((strcmp(opt.method,'AMPhadamardSeeded') == 1 ) || (strcmp(opt.method,'AMPhadamardSeededTranspose') == 1) )
    M = opt.M; N = opt.N;
    if (opt.remove_mean > 0);
        
        vec = zeros(1, opt.M);
        Gmean = zeros(opt.numBlockL, opt.N);
        Ymean = zeros(1, opt.M);
        for l = 1 : opt.numBlockL
            vec(sum(opt.Mblock(1 : (l - 1) ) ) + 1 : sum(opt.Mblock(1 : l) ) ) = 1;
            Gmean(l, :) = 1 / opt.Mblock(l) * MultSeededHadamardTranspose2(vec, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError)';
            Ymean(sum(opt.Mblock(1 : (l - 1) ) ) + 1 : sum(opt.Mblock(1 : l) ) ) = mean(Y(sum(opt.Mblock(1 : (l - 1) ) ) + 1 : sum(opt.Mblock(1 : l) ) ) );
            vec = zeros(1, opt.M);
        end
        
        Gmean2 = Gmean.^2;
        Yeff = Y - Ymean;
    end
end