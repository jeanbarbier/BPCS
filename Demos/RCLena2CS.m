% addpath(genpath('./BPCS'));
clear all
close all

%% Parameters
% size of the patches in which the image is decomposed (power of 2, N >= 64 becomes very demanding in memory)
N = 256.^2;
% do you want to use hadamard matrix (way more fast) or gaussian one?
hadamard = 1;
% seeded or homogeneous matrix (if seeded = 0)? for both there are random and hadamard matrices
seeded = 1;
% measurement rate
alpha = 0.15;
% for seeded matrices, do you want to observe the MSE by block during the reconstruction process every MSEbyBlock iteration? (no = -1)
MSEbyBlock = 10;
% variance of the small additive gaussian noise z
var_1 = 1e-6;
% density of the big noise
rho_corrupt = 0.1;
% seeded decoding matrix parameters (for small alpha)
if (hadamard == 1)
    if (seeded == 1)
        JJ = 0.2; % Correlation parameter between the blocks
        numBlockC = 2^2; % number of blocks for the columns of the seeding matrix (POWER OF 2)
        numBlockL = numBlockC + 1; % number of blocks for the rows of the seeding matrix
        w = 1; % number of sub-diagonal blocks
        alpha1 = rho_corrupt * 2.5;
        alpha2 = (numBlockL * alpha - alpha1) / (numBlockL - 1);
    else
        alpha2 = 0; alpha1 = alpha; JJ = []; numBlockC = 1; numBlockL = 1; w = 0;
    end
else
    if (seeded == 1)
        numBlockC = 10; % number of blocks for the columns of the seeding matrix
        JJ = .02; % Correlation parameter between the blocks
        w = 4; % number of sub-diagonal blocks
        alpha1 = .22; alpha2 = alpha - alpha1 ./ numBlockC;
    end
end

% do you want to use l1magic pack (that must be added to the path) or the inluded l1 AMP solver?
% default is zero (AMP l1 solver), be sure to download l1_magic and put it in the path if you use it: http://users.ece.gatech.edu/~justin/l1magic/
l1_magic = 0; l1_AMP = 1 - l1_magic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation
M = ceil(N ./ (1 - alpha) );
if (hadamard == 1);
    if (seeded == 0)
        NN = 2^(ceil(log(M) / log(2) ) );  % number of variables extended to the first next power of two
        Nblock = NN;
        MblockF = ceil(alpha * Nblock); % number of measurements
        MblockA = ceil((Nblock - MblockF) / numBlockL);
    else
        NN = 2^(ceil(log(M) / log(2) ) ); % number of variables extended to the first next power of two
        Nblock = 2^(ceil(log(NN / numBlockC) / log(2) ) ); % number of variables per block in the seeded matrix (POWER OF 2)
        MblockF(1) = ceil(alpha1 * Nblock); % number of lines of seed block
        MblockF(2 : numBlockL) = ceil(alpha2 * Nblock); % number of variables of bulk blocks
        
        Jf = createSeededHadamardJ(MblockF, Nblock, numBlockL, numBlockC, JJ, w);
        %         MblockA(1) = floor((Nblock - sum(MblockF(:) .* (Jf(:, 1) ~= 0) ) ) / (w + 2) ); % number of lines of seed block
        %         MblockA(2) = floor((Nblock - sum(MblockF(:) .* (Jf(:, 2) ~= 0) ) ) / (w + 3) ); % number of lines of seed block
        %         MblockA(1 : numBlockL) = ceil((Nblock - sum(MblockF(:) .* (Jf(:, 3) ~= 0) ) ) / (w + 3) );
        MblockA(1 : numBlockL) = ceil(0.1 * Nblock);
        
        if (max(MblockA(1) + MblockF(1), MblockA(2) + MblockF(2) ) > Nblock); error('too big A'); end
        
        Ja = createSeededHadamardJ(MblockA, Nblock, numBlockL, numBlockC, JJ, w);
    end
    
    for l = 1 : numBlockL
        vecF = [1 : MblockF(l) ]; vecF = vecF(randperm(MblockF(l) ) );
        vecA = [1 : MblockA(l) ]; vecA = vecA(randperm(MblockA(l) ) );
        noBlockErrorF{l} = vecF(1 : ceil(MblockF(l) / 2) ); % which line are changed of sign (to avoid problems with the reconstruction of the first components of the blocks of the signal)
        noBlockErrorA{l} = vecA(1 : ceil(MblockA(l) / 2) ); % which line are changed of sign (to avoid problems with the reconstruction of the first components of the blocks of the signal)
    end
    
    for c = 1 : numBlockC
        rpNblock = randperm(Nblock);
        rpCf = rpNblock(1 : sum(MblockF .* (Jf(:, c) ~= 0)' ) );
        rpCa = rpNblock(sum(MblockF .* (Jf(:, c) ~= 0)' ) + 1 : end);
        
        u = 1; v = 1;
        for l = 1 : numBlockL
            
            if (Jf(l, c) ~= 0);
                rpermF{l, c} = rpCf(u : u + MblockF(l) - 1); u = u + MblockF(l); % random permutations of the Hadamard modes
                z = 1;
                while (sum(rpermF{l, c}(1 : MblockF(l) ) == 1) > 0); rpermF{l, c} = (rpermF{l, c} ~= 1) .* rpermF{l, c} + (rpermF{l, c} == 1) .* rpermF{l, c}(z); z = z + 1; end % avoid the first (only ones) mode and to have two times the same mode
            end
            
            if (Ja(l, c) ~= 0);
                rpermA{l, c} = rpCa(v : v + MblockA(l) - 1); v = v + MblockA(l); % random permutations of the Hadamard modes
                z = 1;
                while (sum(rpermA{l, c}(1 : MblockA(l) ) == 1) > 0); rpermA{l, c} = (rpermA{l, c} ~= 1) .* rpermA{l, c} + (rpermA{l, c} == 1) .* rpermA{l, c}(z); z = z + 1; end % avoid the first (only ones) mode and to have two times the same mode
            end
        end
    end
end
Lena = imread('Lena256.png');
I_L1 = zeros(size(Lena) ); I_CS = zeros(size(Lena) ); I_id = zeros(size(Lena) );

% Generate the measurement matrix and the coding one
if (hadamard == 0)
    if (seeded == 1)
        [F, Jmat, Mvec, Nvec, alpha] = Create_MAT1D_gauss_Wplus(M, numBlockC, JJ, alpha1, alpha2, w); F = full(F)'; % Gaussian seeded matrix
    else
        F = randn(M,M - N); % Gaussian full matrix
    end
else
%     [F, Jf] = createSeededHadamardMat(MblockF, Nblock, numBlockL, numBlockC, JJ, rpermF, w, noBlockErrorF); % Hadamard matrix;
%     [A, Ja] = createSeededHadamardMat(MblockA, Nblock, numBlockL, numBlockC, JJ, rpermA, w, noBlockErrorA); % Hadamard matrix;
    %     imshow(F)
    %     figure; imshow(A)
    %     max(max(F*A'))
    %     pause
end

disp('Lena decoding')
disp('true gamma')
disp(1 / (1 - sum(MblockF) / NN) )
disp('true alpha')
disp(sum(MblockF) / NN)

for ll = 1 : max(size(Lena) ) / sqrt(N)
    for cc = 1 : max(size(Lena) ) / sqrt(N)
        
        %% Generate the signal and the perfect measure
        S = double(Lena((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1 : cc .* sqrt(N) ) ) ./ 256; % patch of Lena
        
        % haar transform
        Sh = haar_transform(sqrt(N), 5, S);
               
        % selection of the bigger haar coefficients (in absolute value)
        [SsortBig, IndSsortBig] = sort(reshape(Sh, 1, N),'descend');
        [SsortSmall, IndSsortSmall] = sort(reshape(Sh, 1, N),'ascend');
        IndSort = [IndSsortBig(1 : ceil(sum(MblockA) / 2) ), IndSsortSmall(1 : floor(sum(MblockA) / 2) )];
        Sbig = SsortBig(1 : ceil(sum(MblockA) / 2) );
        Ssmall = SsortSmall(1 : floor(sum(MblockA) / 2) );
        SS = zeros(1,N);
        SS(IndSort) = [Sbig, Ssmall];
        SSS = SS(SS ~= 0);
        IndSSS = find(SS ~= 0);

        rp1 = [1:(sum(MblockA) )];
        rp1 = randperm(sum(MblockA) );
        S_mixed = SSS(rp1); % mixing of the components
        if (hadamard == 0); A = null(F'); end % matrix null space
        %         comp = 0; while (min(size(A) ) > max(size(S_mixed) ) ); S_mixed = [S_mixed; 0]; comp = comp + 1; end
        if (hadamard == 0); Y_perfect = A * S_mixed; else Y_perfect = MultSeededHadamardTranspose2(S_mixed, Ja, numBlockL, numBlockC, MblockA, Nblock, rpermA, noBlockErrorA); end
        
        max(max(MultSeededHadamard2(Y_perfect, Jf, numBlockL, numBlockC, MblockF, Nblock, rpermF, noBlockErrorF)))
%         pause
%         Y_perfect2 = A' * S_mixed';
%         max(max(abs(Y_perfect - Y_perfect2) ) )
%         pause
        MM = max(size(Y_perfect) );
        
        %% Error due to noisy channel transmission
        k = floor(rho_corrupt .* MM);
        rp2 = randperm(MM); error = [[randn(k, 1); zeros(MM - k, 1)] ]; error_ = error(rp2); error = error_;
        small_error = randn(MM, 1) .* sqrt(var_1);
        Y_err = Y_perfect + small_error + error;
        if (hadamard == 0); Y_null = F' * Y_err; else Y_null =  MultSeededHadamard2(Y_err, Jf, numBlockL, numBlockC, MblockF, Nblock, rpermF, noBlockErrorF); end % measure
        
%                 Y_null2 = F * Y_err;
%                 max(max(abs(F*A')))
%                 max(abs(MultSeededHadamard2(Y_perfect, Jf, numBlockL, numBlockC, MblockF, Nblock, rpermF, noBlockErrorF) ) )
%                 max(abs(Y_null2 - Y_null))
%                 pause
        
        %         %% Ideal coding and decoding
        %         Y_id = Y_perfect + small_error;
        %         X_id = A \ Y_id; X_id_(rp) = X_id(1:max(size(X_id) ) - comp, 1); X_id = X_id_;
        %
        %         %% Naive decoding (with no error correction)
        %         X_naive = A \ Y_err; X_naive_(rp) = X_naive(1:max(size(X_naive) ) - comp, 1); X_naive = X_naive_;
        
        %% BPCS coding and decoding by reconstruction of the sparse error vector
        % algorithm properties
        My = CSBP_Solver_Opt();
        My.nb_iter = 1000;
        My.print = 10;
        My.conv = var_1 * 5;
        My.signal_rho = rho_corrupt;
        My.signal = error + small_error;
        My.dump_mes = 0.5;
        My.MSEbyBlock = MSEbyBlock;
        My.noBlockError = noBlockErrorF;
        if (hadamard == 1);
            My.method = 'AMPhadamardSeeded';
            My.N = NN;
            My.M = sum(MblockF);
            My.J = Jf;
            My.numBlockL = numBlockL;
            My.numBlockC = numBlockC;
            My.Mblock = MblockF;
            My.Nblock = Nblock;
            My.rp = rpermF;
        else My.method = 'AMP';
        end
        % 2Gauss
        My.m_2_gauss = 0;
        My.var_1_gauss = var_1;
        My.var_2_gauss = 1;
        
        %% AMP decoding with 2Gauss prior
        My.prior = '2Gauss';
        if (hadamard == 0); [results2G, n_and_e] = CSBP_Solver(Y_null', F', My); else [results2G, n_and_e] = CSBP_Solver(Y_null', [], My); end;
        Y_CS = Y_err - results2G.av_mess';
        pause
        
        if (hadamard == 0);
            X_CS(rp1) = A \ Y_CS;
        else
            % algorithm properties
            My = CSBP_Solver_Opt();
            My.nb_iter = 100;
            My.print = 10;
            My.conv = 1e-10;
            My.signal_rho = 1;
            My.signal = S_mixed;
            My.dump_mes = 0.;
            My.MSEbyBlock = MSEbyBlock;
            My.noBlockError = noBlockErrorA;
            My.method = 'AMPhadamardSeededTranspose';
            My.N = sum(MblockA);
            My.M = NN;
            My.J = Jf;
            My.numBlockL = numBlockL;
            My.numBlockC = numBlockC;
            My.Mblock = MblockA;
            My.Nblock = Nblock;
            My.rp = rpermA;
            % SparseGauss
            My.prior = 'SparseGauss';
            My.m_gauss = mean(SSS);
            My.var_gauss = var(SSS);
            [resultsCS, n_and_e] = CSBP_Solver(Y_CS, [], My);
        end
        X_CS(rp1) = resultsCS.av_mess;
        S_CS = zeros(1, N);
        S_CS(IndSSS) = X_CS;
        SS_CS = haar_inverse(sqrt(N), 5, reshape(S_CS, sqrt(N), sqrt(N) ) );
        imshow(SS_CS)
        pause
        
        
        %% L1
        if (l1_AMP == 1)
            My.prior = 'L1';
            [resultsL1 ,n_and_e] = CSBP_Solver(Y_null', F', My);
            Y_L1 = Y_err - resultsL1.av_mess';
            X_L1 = A \ Y_L1; X_L1_(rp) = X_L1(1 : max(size(X_L1) ) - comp, 1); X_L1 = X_L1_;
        else
            % Dantzig selector
            err_L1 = l1dantzig_pd(F * Y_null, F', [], Y_null, 3e-3, 1e-2, 1000);
            % regression step
            err_index = abs(err_L1) > sqrt(var_1 .* M); F_t = F';
            F_index_t = F_t(:,err_index);
            err_L1 = zeros(size(Y_err) ); err_L1(err_index) = F_index_t \ Y_null;
            Y_L1 = Y_err - err_L1;
            X_L1 = A \ Y_L1; X_L1_(rp) = X_L1(1 : max(size(X_L1) ) - comp, 1);  X_L1 = X_L1_;
        end
        
        %% Results
        S_unmixed(rp) = S(1 : max(size(S_mixed) ) - comp, 1); S = S_unmixed;
        rho_id_CS = norm(X_CS - S, 2) ./ norm(X_id - S, 2);
        rho_id_L1 = norm(X_L1 - S, 2) ./ norm(X_id - S, 2);
        mse_id = mean((X_id - S).^2);
        mse_CS = mean((X_CS - S).^2);
        mse_L1 = mean((X_L1 - S).^2);
        
        results_(1, ll, cc) = {reshape(X_id, sqrt(N), sqrt(N) )};
        results_(2, ll, cc) = {reshape(X_CS, sqrt(N), sqrt(N) )};
        results_(3, ll, cc) = {reshape(X_L1, sqrt(N), sqrt(N) )};
        results_(4, ll, cc) = {[mse_id, mse_CS, mse_L1, rho_id_CS, rho_id_L1]};
        
        I_CS((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1: cc .* sqrt(N) ) = reshape(X_CS, sqrt(N), sqrt(N) );
        I_L1((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1: cc .* sqrt(N) ) = reshape(X_L1, sqrt(N), sqrt(N) );
        I_id((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1: cc .* sqrt(N) ) = reshape(X_id, sqrt(N), sqrt(N) );
        I_naive((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1: cc .* sqrt(N) ) = reshape(X_naive, sqrt(N), sqrt(N) );
        
    end
end

%% Plots
clf
subplot(3, 2, 1); imshow(Lena); xlabel('Original Image');
subplot(3, 2, 2); imshow(I_id); xlabel('Ideal');
subplot(3, 2, 3); imshow(I_CS); xlabel('AMP 2Gauss');
subplot(3, 2, 4); imshow(I_L1); xlabel('L1');
subplot(3, 2, 5); imshow(I_naive); xlabel('No error correction');

%% Save
save('Lena_decoding','I_id','I_CS','I_L1','I_naive','results_');