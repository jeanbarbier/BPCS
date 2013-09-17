close all
new = 1; % new = 1 for a new CS (Signal + matrix) instance
if (new == 1); clear; new = 1; end
realMat = 0; % realMat = 1 if you want to actually create and see the seeded matrix (FOR NOT TOO BIG N)

% properties
choice = 'SparseGauss'; % Choice of the signal and the associated prior
rho = 0.1; % density signal
N = 2^15; % Size of signal POWER OF 2;
JJ = 0.3; % Correlation parameter of the upper diagonal blocks
numBlockC = 2^4; % number of blocks for the columns of the seeding matrix (POWER OF 2)
numBlockL = numBlockC + 1; % number of blocks for the rows of the seeding matrix
Nblock = N / numBlockC; % number of variables per block in the seeding matrix (POWER OF 2)
Mblock(1) = ceil(Nblock / 4); % number of measurements of the first (seed) block
Mblock(2 : numBlockL) = ceil(Nblock / 7); % number of measurements of the next ones
w = 3; % number of sub-diagonal blocks

alpha1 = Mblock(1) / Nblock % alpha seed
alpha2 = Mblock(2) / Nblock % alpha bulk
alpha = sum(Mblock) / N % alpha total

G = [];
if (new == 1)
    
    % signal
    switch choice
        case 'SparseGauss'
            average = 0; var = 1; S = S_SparseGauss(N,rho,average,var);
        case 'SparseExponential'
            exposent = 1; S = S_SparseExponential(N,rho,exposent);
        case 'SparseBinary'
            S = S_SparseBinary(size_S,rho);
        case '2Gauss'
            m_1 = 0; m_2 = 1; var_1 = 1e-8; var_2 = 1; S = S_2Gauss(N,rho,m_1,m_2,var_1,var_2);
        case 'SparseConstant'
            down = -1; up = 1;  S = S_SparseConstant(N,rho,down,up);
        case 'L1'
            average = 0; var = 1; S = S_SparseGauss(N,rho,average,var);
        case 'Laplace'
            S = (-log(rand(N,1)).*(2.*(rand(N,1)>0.5)-1));
        otherwise
            disp('Error in choice')
    end
    
    for l = 1 : numBlockL
        vec = [1 : Mblock(l) ];
        vec = vec(randperm(Mblock(l) ) );
        noBlockError{l} = vec(1 : ceil(Mblock(l) / 4) ); % which line are changed of sign (to avoid problems with the reconstruction of the first components of the blocks of the signal)
        for c = 1 : numBlockC
            rp{l, c} = randperm(Nblock); % random permutations of the Hadamard modes
            if (sum(rp{l, c}(1 : Mblock(l) ) ~= 1) > 0);
                rp{l, c} = (rp{l, c} ~= 1) .* rp{l, c} + (rp{l, c} == 1) .* rp{l, c}(end); % avoid the first (only ones) mode and to have two times the same mode
            end
        end
    end
    
    if (realMat == 1)
        % matrix creation
        [G, J] = createSeededHadamardMat(Mblock, Nblock, numBlockL, numBlockC, JJ, rp, w, noBlockError);
        figure; imshow(G); drawnow;
        
        % measure
        Y = G * S';
        
    else
        
        J = createSeededHadamardJ(Mblock, Nblock, numBlockL, numBlockC, JJ, w); % Create the J (variance by block) matrix
        
        % measure without the matrix in memory (allows way bigger N)
        Y = MultSeededHadamard2(S, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError);
        
    end
    
end

% tests des differents opérateurs
% a = randn(N, 1);
% u = G * a;
% uu = MultSeededHadamard2(a, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError);
% figure; plot(abs(u - uu) );
%
% b = randn(1, sum(Mblock) );
% v = b * G;
% vv = MultSeededHadamardTranspose2(b, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)';
% figure; plot(abs(v - vv) );
%
% a = randn(N, 1);
% u = G.^2 * a;
% uu = MultSeededHadamardSquarred2(a, J, numBlockL, numBlockC, Mblock, Nblock);
% figure; plot(abs(u - uu) );
%
% b = randn(1, sum(Mblock) );
% v = b * G.^2;
% vv = MultSeededHadamardTransposeSquarred2(b, J, numBlockL, numBlockC, Mblock, Nblock)';
% figure; plot(abs(v - vv) );

% algo properties
My = CSBP_Solver_Opt();
My.MSEbyBlock = 20; % Show the MSE by blocks of the signal, -1 for no use
My.nb_iter = 5000;
My.print = 1;
My.conv = 10^(-10);
My.learn = 0;
My.signal_rho = -1; % (-1 -> measure_rate ./ 10)
if (My.learn == 0); My.signal_rho = rho; end;
My.option_noise = 0;
My.var_noise = 0;
My.signal = S;
My.dump_learn = 0.5;
My.dump_mes = 0.5;
My.prior = choice;
My.remove_mean = 1; % must remain 0 for the moment
My.method = 'AMPhadamardSeeded';
% My.method = 'AMP';
My.M = sum(Mblock);
My.N = N;
My.J = J;
My.numBlockL = numBlockL;
My.numBlockC = numBlockC;
My.Mblock= Mblock;
My.Nblock = Nblock;
My.rp = rp;
My.noBlockError = noBlockError;
% SparseGauss
if ((My.learn == 0) && (strcmp(choice,'SparseGauss') ) )
    My.m_gauss = average;
    My.var_gauss = var;
else
    My.m_gauss = 0;
    My.var_gauss = 1;
end
% 2Gauss
if ((My.learn == 0) && (strcmp(choice,'2Gauss') ) )
    My.m_1_gauss = m_1;
    My.m_2_gauss = m_2;
    My.var_1_gauss = var_1;
    My.var_2_gauss = var_2;
else
    My.m_1_gauss = 1;
    My.m_2_gauss = 1;
    My.var_1_gauss = 1e-5;
    My.var_2_gauss = 1;
end
% SparseExponential
if ((My.learn == 0) && (strcmp(choice,'SparseExponential') ) )
    My.expo = exposent;
else
    My.expo = 1;
end
% SparseConstant
My.c_down = -1;
My.c_up = 1;
% L1
My.min = -10;
My.max = 10;
% Laplace
My.beta = 1;

% algo
figure;
[results,n_and_e] = CSBP_Solver(Y,G,My);