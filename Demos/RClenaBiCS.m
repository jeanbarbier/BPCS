addpath(genpath('~/Documents/taff/CS/RC_Demos3') ); clear all; close all;

%% Parameters
% size of Lena
N = 64^2;
% fraction of wavelet coefficients that we keep/aim to reconstruct
alphaWavelet = 0.06;
% coding rate (power of 2)
M = N / 2^2;
% measurement rate
alphaCS = 0.04;
% do you want to use hadamard matrix (way more fast) or gaussian one?
hadamard = 1;
% seeded or homogeneous matrix (if seeded = 0)? for both there are random and hadamard matrices
seeded = 1;
% for seeded matrices, do you want to observe the MSE by block during the reconstruction process every MSEbyBlock iteration? (no = -1)
MSEbyBlock = 10;
% variance of the small additive gaussian noise z
var_1 = 1e-10;
% density of the big noise
rho_corrupt = 0.01;
% number of haar transform of the original picture (<= 8 for a N = 256^2 picture)
levelHaar = 6;
% Do you want to generate the seeded coding and decoding matrices even if not required? (for smaller N)
realMat = 1;
% if realMat = 1; do you want to visualize the matrices?
showMat = 0;
% seeded decoding matrix parameters (for small alpha)
if (hadamard == 1)
    if (seeded == 1)
        % Correlation parameter between the blocks
        JJ = 0.2;
        % number of blocks for the columns of the seeding matrix (POWER OF 2)
        numBlockC = 16;
        % number of blocks for the rows of the seeding matrix
        numBlockL = numBlockC + 1;
        % number of blocks for the rows of the coding matrix
        numBlockL_A = 8 * numBlockC;
        % number of sub-diagonal blocks
        w = 1;
        % alpha 1st block (seed)
        alphaCS1 = rho_corrupt * 2.5;
        % alpha bulk blocks
        alphaCS2 = (numBlockL * alphaCS - alphaCS1) / (numBlockL - 1);
    else
        alphaCS2 = 0; alphaCS1 = alphaCS; JJ = []; numBlockC = 1; numBlockL = 1; w = 0;
    end
else
    if (seeded == 1)
        % number of blocks for the columns of the seeding matrix
        numBlockC = 10;
        % Correlation parameter between the blocks
        JJ = .02;
        % number of sub-diagonal blocks
        w = 4;
        % alpha 1st block (seed)
        alphaCS1 = .22;
        % alpha bulk blocks
        alphaCS2 = alphaCS - alphaCS1 ./ numBlockC;
    end
end

% check on the size not to saturate memory
if ((realMat == 1) && (sqrt(N) > 64) ); error('Too big N to create the measurement and coding matrices'); end

% do you want to use l1magic pack (that must be added to the path) or the inluded l1 AMP solver?
% default is zero (AMP l1 solver), be sure to download l1_magic and put it in the path if you use it: http://users.ece.gatech.edu/~justin/l1magic/
l1_magic = 0; l1_AMP = 1 - l1_magic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation
if (hadamard == 1);
    if (seeded == 0)
        Nblock = M;
        % number of measurements of the F matrix
        MblockF = ceil(alphaCS * Nblock);
        % number of measurements of the A matrix
        MblockA = ceil((Nblock - MblockF) / numBlockL);
    else
        % number of variables per block in the seeded matrix (POWER OF 2)
        Nblock = 2^(ceil(log(M / numBlockC) / log(2) ) );
        % number of lines of seed block
        MblockF(1) = floor(alphaCS1 * Nblock);
        % number of variables of bulk blocks
        MblockF(2 : numBlockL) = floor(alphaCS2 * Nblock);
        % number of variables of bulk blocks for the coding matrix
        MblockA(1 : numBlockL_A) = ceil(N / numBlockL_A);
        % creation of the matrices containing the amplitudes by block of the measurement matrix F
        Jf = createSeededHadamardJ(numBlockL, numBlockC, JJ, w);
        % creation of the matrices containing the amplitudes by block of the measurement matrix A
        Ja = randn(numBlockL_A, numBlockC);
    end
    
    for l = 1 : numBlockL
        vecF = [1 : MblockF(l) ]; vecF = vecF(randperm(MblockF(l) ) );
        % which line are changed of sign (to avoid problems with the reconstruction of the first components of the blocks of the signal)
        noBlockErrorF{l} = vecF(1 : ceil(MblockF(l) / 2) );
    end
    
    for l = 1 : numBlockL_A
        vecA = [1 : MblockA(l) ]; vecA = vecA(randperm(MblockA(l) ) );
        % which line are changed of sign (to avoid problems with the reconstruction of the first components of the blocks of the signal)
        noBlockErrorA{l} = vecA(1 : ceil(MblockA(l) / 2) );
    end
    
    % computation of the fraction of Hadamard modes used for the seeded matrix
    borne = 0;
    for c = 1 : numBlockC
        borne = max(borne, sum(MblockF' .* (Jf(:, c) ~= 0) ) );
    end
    disp('fraction of Hadamard modes used for the seeded matrix'); disp(borne / Nblock);
    
    for c = 1 : numBlockC
        fullPerm = randperm(Nblock);
        i1 = 1; while (i1 <= ceil(Nblock / 2) ); fullPerm = randperm(Nblock); i1 = find(fullPerm == 1); end
        % No first mode
        fullPerm(i1) = fullPerm(end);
        
        permF = fullPerm(1 : borne);
        permA = fullPerm(borne + 1 : end - 1);
        if (max(size(permA) ) < max(MblockA) ); error('You must increase numBlockL_A to avoid repetitions of Hadamard modes in a block'); end
        
        rpCf = permF(randperm(max(size(permF) ) ) );
        rpCa = permA(randperm(max(size(permA) ) ) );
        
        u = 1;
        for l = 1 : numBlockL
            if (Jf(l, c) ~= 0); rpermF{l, c} = rpCf(u : u + MblockF(l) - 1); u = u + MblockF(l);
            else rpermF{l, c} = []; end
        end
        
        v = 1;
        for l = 1 : numBlockL_A
            if (Ja(l, c) ~= 0);
                if (v >= max(size(rpCa) ) - MblockA(1) ); v = 1; rpCa = permA(randperm(max(size(permA) ) ) ); end
                rpermA{l, c} = rpCa(v : v + MblockA(l) - 1); v = v + MblockA(l);
            else rpermA{l, c} = []; end
        end
    end
end

Lena = imread('Lena256.png');
I_L1 = zeros(size(Lena) ); I_CS = zeros(size(Lena) ); I_id = zeros(size(Lena) );

% Generate the measurement matrix and the coding one
if (hadamard == 0)
    if (seeded == 1)
        % Gaussian seeded matrix
        [F, Jmat, Mvec, Nvec, alpha] = Create_MAT1D_gauss_Wplus(M, numBlockC, JJ, alphaCS1, alphaCS2, w); F = full(F)';
    else
        % Gaussian full matrix
        F = randn(M,M - N);
    end
else
    if (realMat == 1)
        % Hadamard matrix;
        [A] = createCodingHadamardMat(MblockA, Nblock, numBlockL_A, numBlockC, rpermA, Ja, noBlockErrorA);
        [F] = createCodingHadamardMat(MblockF, Nblock, numBlockL, numBlockC, rpermF, Jf, noBlockErrorF);
        if (showMat == 1); figure; imshow(F); figure; imshow(A); pause; end
    end
end

disp('Lena decoding');
disp('true gamma'); disp(M / N);
disp('true alpha'); disp(sum(MblockF) / M);

for ll = 1 : max(size(Lena) ) / sqrt(N)
    for cc = 1 : max(size(Lena) ) / sqrt(N)
        
        %% Generate the signal and the perfect measure
        % patch of Lena
        S = double(Lena((ll - 1) * sqrt(N) + 1: ll * sqrt(N), (cc - 1) * sqrt(N) + 1 : cc * sqrt(N) ) ) ./ 256;
        
        % haar transform
        Sh = haar_transform(sqrt(N), levelHaar, S);
        
        % selection of the bigger haar coefficients (in absolute value)
        [SsortBig, IndSsortBig] = sort(reshape(Sh, 1, N),'descend');
        [SsortSmall, IndSsortSmall] = sort(reshape(Sh, 1, N),'ascend');
        IndSort = [IndSsortBig(1 : ceil(alphaWavelet * N / 2) ), IndSsortSmall(1 : ceil(alphaWavelet * N / 2) )];
        Sbig = SsortBig(1 : ceil(alphaWavelet * N / 2) );
        Ssmall = SsortSmall(1 : ceil(alphaWavelet * N / 2) );
        SS = zeros(1, N); SS(IndSort) = [Sbig, Ssmall];
        
        % mixing of the components
        rp1 = randperm(N); S_mixed = SS(rp1);
        S_mixed = [S_mixed, zeros(1, sum(MblockA) - N) ];
        disp('number of concatenated zeros to the true signal'); disp(sum(MblockA) - N);
        
        % matrix null space
        if (hadamard == 0); A = null(F'); Y_perfect = A * S_mixed;
        else Y_perfect = MultSeededHadamardTransposeA(S_mixed, Ja, numBlockL_A, numBlockC, MblockA, Nblock, rpermA, noBlockErrorA); end
        M = max(size(Y_perfect) );
        
%         A = randn(size(A));
%         Y_perfect = A' * S_mixed';
        
        % test of the operators
        if (realMat == 1)
            disp('***** beginning test of the operators, you should observe ~ 0 *****');
            max(Y_perfect - A' * S_mixed')
            SM = randn(M, 1);
            SN = randn(1, sum(MblockA) );
            test1 = MultSeededHadamardTransposeA(SN, Ja, numBlockL_A, numBlockC, MblockA, Nblock, rpermA, noBlockErrorA);
            test_1 = SN * A;
            max(abs(test1' - test_1) )
            test2 = MultSeededHadamardTransposeSquarredA(SN, Ja, numBlockL_A, numBlockC, MblockA, Nblock);
            test_2 = SN * A.^2;
            max(abs(test2' - test_2) )
            test3 = MultSeededHadamardA(SM, Ja, numBlockL_A, numBlockC, MblockA, Nblock, rpermA, noBlockErrorA);
            test_3 = A * SM;
            max(abs(test3 - test_3) )
            test4 = MultSeededHadamardSquarredA(SM, Ja, numBlockL_A, numBlockC, MblockA, Nblock);
            test_4 = A.^2 * SM;
            max(abs(test4 - test_4) )
        end
        
        
        %% Error due to noisy channel transmission
        k = floor(rho_corrupt .* M);
        rp2 = randperm(M); error = [[randn(k, 1); zeros(M - k, 1) ] ]; error = error(rp2);
        small_error = randn(M, 1) * sqrt(var_1);
        Y_err = Y_perfect + small_error + error;
        
        % measure
        if (hadamard == 0); Y_null = F' * Y_err; else Y_null =  MultSeededHadamard2(Y_err, Jf, numBlockL, numBlockC, MblockF, Nblock, rpermF, noBlockErrorF); end
        
        % additional tests of the operators
        if (realMat == 1)
            Y_null2 = F * Y_err;
            max(abs(Y_null2 - Y_null) )
            max(max(abs(F * A') ) )
            max(abs(MultSeededHadamard2(Y_perfect, Jf, numBlockL, numBlockC, MblockF, Nblock, rpermF, noBlockErrorF) ) )
            disp('***** end of the test of the operators *****');
        end
        
        %% Ideal decoding and decoding
        disp('ideal decoding')
        Y_id = Y_perfect + small_error;
        
        %         As = randn(size(A) );
        %         Ar = (As >= 0) - (As < 0);
        %
        %         Ss = randn(size(S_mixed) );
        %         Sr = (Ss >= 0) - (Ss < 0);
        %         S_mixedd = S_SparseGauss(max(size(S_mixed)),sum(S_mixed ~= 0) / max(size(S_mixed) ), mean(S_mixed) / (sum(S_mixed ~= 0) / max(size(S_mixed) ) ), var(S_mixed) / (sum(S_mixed ~= 0) / max(size(S_mixed) ))^2);
        %         S_mixedd = Sr .* S_mixed;
        
        %         imshow(A); drawnow; pause
        
        %         [lin, col] = size(A);
        %         for c = 1 : col
        %             ranC(c) = randn;
        %             A(:, c) = ranC(c) .* A(:, c);
        %         end
        %         for l = 1 : lin
        %             ranL(1, l) = randn;
        %             A(l, :) = ranL(l) .* A(l, :);
        %         end
        
        %         srpTryL = randperm(lin);
        %         rpTryC = randperm(col);
        %         A = A(rpTryL, rpTryC);
        %         rpermA = rpermA{rpTry, :};
        %                 for l = 1 : lin
        %                     for c = 1 : col
        %                         ran1 = randn;
        %                         ran2 = randn;
        %                         ranLCind(l, c) = ran1 * ran2;
        %                         ranLCcor(l, c) = ranC(c) * ranL(l);
        %                         oneLCind(l, c) = ran1 * ran2 > 0;
        %                         oneLCcor(l, c) = ranC(c) * ranL(l) > 0;
        %                     end
        %                 end
        
        
        
        %         for c = 1 : col
        %             for l = 1 : lin
        %                 if (ranC(c) * ranL(l) > 0); A(l, c) = -A(l, c); end
        %             end
        %         end
        %         imshow(A); drawnow; pause
        
        
        %         Y_id =  A' * S_mixed';
        
        
        if (hadamard == 0);
            X_id(rp1) = A \ Y_id;
        else
            % algorithm properties
            My = CSBP_Solver_Opt();
            My.nb_iter = 20000;
            My.print = 1;
            My.conv = 1e-8;
            My.signal_rho = sum(S_mixed ~= 0) / max(size(S_mixed) );
            My.signal = S_mixed;
            My.dump_mes = 0.5;
            My.MSEbyBlock = MSEbyBlock;
            My.noBlockError = noBlockErrorA;
            My.method = 'AMPhadamardSeededTransposeA';
            if (realMat == 1); My.method = 'AMP'; end
            My.remove_mean = 0;
            My.learn = 0;
            My.dump_learn = 0.5;
            My.M = N;
            My.N = M;
            My.J = Ja;
            My.numBlockL = numBlockL_A;
            My.numBlockC = numBlockC;
            My.Mblock = MblockA;
            My.Nblock = Nblock;
            My.rp = rpermA;
            % SparseGauss
            My.prior = 'SparseGauss';
            My.m_gauss = mean(S_mixed) / My.signal_rho;
            My.var_gauss = var(S_mixed) / My.signal_rho^2;
            if (realMat == 1); [resultsCS, n_and_e] = CSBP_Solver(Y_id, A', My);
            else [resultsCS, n_and_e] = CSBP_Solver(Y_id, [], My); end
            X_id(rp1) = resultsCS.av_mess(1 : N);
        end
        %         S_id = zeros(1, N); S_id(IndSSS) = X_id;
        %         X_id = haar_inverse(sqrt(N), levelHaar, reshape(S_id, sqrt(N), sqrt(N) ) );
        X_id = haar_inverse(sqrt(N), levelHaar, reshape(X_id, sqrt(N), sqrt(N) ) );
        pause
        
        
        %% Naive decoding (with no error correction)
        %         disp('naive decoding')
        %         if (hadamard == 0);
        %             X_naive(rp1) = A \ Y_err;
        %         else
        %             [resultsCS, n_and_e] = CSBP_Solver(Y_err, [], My);
        %             X_naive(rp1) = resultsCS.av_mess;
        %         end
        %         S_naive = zeros(1, N); S_naive(IndSSS) = X_naive;
        %         X_naive = haar_inverse(sqrt(N), levelHaar, reshape(S_naive, sqrt(N), sqrt(N) ) );
        
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
            My.N = M;
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
        disp('AMP estimation of the error')
        My.prior = '2Gauss';
        if (hadamard == 0); [results2G, n_and_e] = CSBP_Solver(Y_null', F', My); else [results2G, n_and_e] = CSBP_Solver(Y_null', [], My); end;
        Y_CS = Y_err - results2G.av_mess';
        
        disp('AMP decoding')
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
            My.dump_mes = 0.9;
            My.MSEbyBlock = MSEbyBlock;
            My.noBlockError = noBlockErrorA;
            My.method = 'AMPhadamardSeededTranspose';
            My.N = sum(MblockA);
            My.M = M;
            My.J = Jf;
            My.numBlockL = numBlockL;
            My.numBlockC = numBlockC;
            My.Mblock = MblockA;
            My.Nblock = Nblock;
            My.rp = rpermA;
            % SparseGauss
            My.prior = 'SparseGauss';
            My.m_gauss = mean(SS);
            My.var_gauss = var(SS);
            [resultsCS, n_and_e] = CSBP_Solver(Y_CS, [], My);
            X_CS(rp1) = resultsCS.av_mess;
        end
        S_CS = zeros(1, N); S_CS(IndSSS) = X_CS;
        X_CS = haar_inverse(sqrt(N), levelHaar, reshape(S_CS, sqrt(N), sqrt(N) ) );
        
        
        %% L1
        if (l1_AMP == 1)
            disp('L1 decoding')
            if (hadamard == 0);
                X_L1(rp1) = A \ Y_CS;
            else
                % algorithm properties
                My.prior = 'L1';
                [resultsCS, n_and_e] = CSBP_Solver(Y_CS, [], My);
                X_L1(rp1) = resultsCS.av_mess;
            end
            S_L1 = zeros(1, N); S_CS(IndSSS) = X_L1;
            X_L1 = haar_inverse(sqrt(N), levelHaar, reshape(S_L1, sqrt(N), sqrt(N) ) );
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
% clf
% subplot(3, 2, 1); imshow(Lena); xlabel('Original Image');
% subplot(3, 2, 2); imshow(I_id); xlabel('Ideal');
% subplot(3, 2, 3); imshow(I_CS); xlabel('AMP 2Gauss');
% subplot(3, 2, 4); imshow(I_L1); xlabel('L1');
% subplot(3, 2, 5); imshow(I_naive); xlabel('No error correction');
%
% %% Save
% save('Lena_decoding','I_id','I_CS','I_L1','I_naive','results_');