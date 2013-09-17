%% Tree path generation

addpath(genpath('~/Documents/taff/CS/RC_Demos3') ); clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General Parameters

% image selection ('Lena' (256^2) or 'Phantom' (any size) )
imageName = 'Lena';
% size of image (power of 2)
N = 64^2;
% fraction of wavelet coefficients that we keep/aim to reconstruct
alphaWavelet = 0.04;
% coding rate (power of 2)
M = N / 2^2;
% measurement rate
alphaCs = 0.1;
% do you want to use hadamard matrix (way more fast) or gaussian one?
hadamard = 1;
% seeded or homogeneous matrix (if seeded = 0)? for both there are random and hadamard matrices
seeded = 1;
% for seeded matrices, do you want to observe the MSE by block during the reconstruction process every MSEbyBlock iteration? (no = -1)
MSEbyBlock = 10;
% variance of the big additive gaussian noise z
varBig = 1e3;
% variance of the small additive gaussian noise z
varSmall = 1e-6;
% density of the big noise
rhoCorrupt = 0.01;
% number of haar transform of the original picture (<= 8 for a N = 256^2 picture)
levelHaar = 6;
% Do you want to generate the seeded coding and decoding matrices even if not required? (for smaller N)
realMat = 1;
% if realMat = 1; do you want to visualize the matrices?
showMat = 1;
% selection of the reconstruction types (0/1)
recIdeal = 1;
recNaive = 1;
recAmp = 1;
recL1 = 1;
% save the results? (0/1)
saveRes = 0;
% name of the file saved
fileName = [imageName, 'Decoding'];
% seeded decoding matrix parameters (for small alpha)
if (hadamard == 1)
    if (seeded == 1)
        % Correlation parameter between the blocks
        JJ = 0.5;
        % number of blocks for the columns of the seeding matrix (POWER OF 2)
        numBlockC = 4;
        % number of blocks for the rows of the seeding matrix
        numBlockL = numBlockC + 1;
        % number of blocks for the rows of the coding matrix
        numBlockLa = 8 * numBlockC;
        % number of sub-diagonal blocks
        w = 2 ;
        % alpha 1st block (seed)
        alphaCs1 = 2 * alphaCs;
        % alpha bulk blocks
        alphaCs2 = (numBlockL * alphaCs - alphaCs1) / (numBlockL - 1);
    else
        alphaCs2 = 0; alphaCs1 = alphaCs; JJ = []; numBlockC = 1; numBlockL = 1; w = 0;
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
        alphaCs1 = .22;
        % alpha bulk blocks
        alphaCs2 = alphaCs - alphaCs1 ./ numBlockC;
    end
end
% general algorithm properties
My = CSBP_Solver_Opt();
My.MSEbyBlock = MSEbyBlock;
% number of maximum iterations
nbIterIdeal = 500;
nbIterNaive = 500;
nbIterAmp = 500;
nbIterL1 = 500;
% priting frequency of results
My.print = 30;
% convergence criterion
My.conv = varSmall * 5;
% dumping of messages
dumpMesA = 0.9;
dumpMesF = 0.9;
% does the coding or decoding matrices have non-zero mean? if yes, remove_mean = 2
My.remove_mean = 0;
% learning ?
My.learn = 0;
% dump of learning
My.dump_learn = 0.5;
% do you want to use l1magic pack (that must be added to the path) or the inluded l1 AMP solver?
% default is zero (AMP l1 solver), be sure to download l1_magic and put it in the path if you use it: http://users.ece.gatech.edu/~justin/l1magic/
l1_magic = 0; l1_AMP = 1 - l1_magic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation and matrices parameters

% check on the size not to saturate memory
if ((realMat == 1) && (sqrt(N) > 64) ); error('Too big N to create the measurement and coding matrices'); end

% matrices parameters
if (hadamard == 1);
    if (seeded == 0)
        Nblock = M;
        % number of measurements of the F matrix
        MblockF = ceil(alphaCs * Nblock);
        % number of measurements of the A matrix
        MblockA = ceil((Nblock - MblockF) / numBlockL);
    else
        % number of variables per block in the seeded matrix (POWER OF 2)
        Nblock = 2^(ceil(log(M / numBlockC) / log(2) ) );
        % number of lines of seed block
        MblockF(1) = floor(alphaCs1 * Nblock);
        % number of variables of bulk blocks
        MblockF(2 : numBlockL) = floor(alphaCs2 * Nblock);
        % number of variables of bulk blocks for the coding matrix
        MblockA(1 : numBlockLa) = ceil(N / numBlockLa);
        % creation of the matrices containing the amplitudes by block of the measurement matrix F
        Jf = createSeededHadamardJ(numBlockL, numBlockC, JJ, w);
        % creation of the matrices containing the amplitudes by block of the measurement matrix A
        Ja = randn(numBlockLa, numBlockC);
    end
    
    % selection of the lines of the matrices that are sign-switched
    for l = 1 : numBlockL
        vecF = [1 : MblockF(l) ]; vecF = vecF(randperm(MblockF(l) ) );
        % which line are changed of sign (to avoid problems with the reconstruction of the first components of the blocks of the signal)
        noBlockErrorF{l} = vecF(1 : ceil(MblockF(l) / 2) );
    end
    
    for l = 1 : numBlockLa
        vecA = [1 : MblockA(l) ]; vecA = vecA(randperm(MblockA(l) ) );
        % which line are changed of sign (to avoid problems with the reconstruction of the first components of the blocks of the signal)
        noBlockErrorA{l} = vecA(1 : ceil(MblockA(l) / 2) );
    end
    
    % computation of the fraction of Hadamard modes used for the seeded matrix
    borne = 0;
    for c = 1 : numBlockC
        borne = max(borne, sum(MblockF' .* (Jf(:, c) ~= 0) ) );
    end
    borneTrue = borne;
    disp('fraction of Hadamard modes used for the seeded matrix'); disp(borneTrue / Nblock);
    
    % selection of the Hadamard modes for each column of the matrices
    for c = 1 : numBlockC
        fullPerm = randperm(Nblock);
        i1 = 1; while (i1 <= ceil(Nblock / 2) ); fullPerm = randperm(Nblock); i1 = find(fullPerm == 1); end
        % No first mode
        fullPerm(i1) = fullPerm(end);
        
        permF = fullPerm(1 : borneTrue);
        permA = fullPerm(borneTrue + 1 : end - 1);
        if (max(size(permA) ) < max(MblockA) ); error('You must increase numBlockL_A to avoid repetitions of Hadamard modes in a block'); end
        
        rpCf = permF(randperm(max(size(permF) ) ) );
        rpCa = permA(randperm(max(size(permA) ) ) );
        
        u = 1;
        for l = 1 : numBlockL
            if (Jf(l, c) ~= 0);
                if (u >= max(size(rpCf) ) - MblockF(1) ); u = 1; rpCf = permF(randperm(max(size(permF) ) ) ); end
                rpermF{l, c} = rpCf(u : u + MblockF(l) - 1); u = u + MblockF(l);
            else rpermF{l, c} = []; end
        end
        
        v = 1;
        for l = 1 : numBlockLa
            if (Ja(l, c) ~= 0);
                if (v >= max(size(rpCa) ) - MblockA(1) ); v = 1; rpCa = permA(randperm(max(size(permA) ) ) ); end
                rpermA{l, c} = rpCa(v : v + MblockA(l) - 1); v = v + MblockA(l);
            else rpermA{l, c} = []; end
        end
    end
end

% image loading
if (strcmp(imageName, 'Lena') ); image = imread('Lena256.png'); else image = phantom(sqrt(N) ); end
I_L1 = zeros(size(image) ); I_CS = zeros(size(image) ); I_id = zeros(size(image) );

% Generate the measurement matrix and the coding one
if (hadamard == 0)
    if (seeded == 1)
        % Gaussian seeded matrix
        [F, Jmat, Mvec, Nvec, alpha] = Create_MAT1D_gauss_Wplus(M, numBlockC, JJ, alphaCs1, alphaCs2, w); F = full(F)';
    else
        % Gaussian full matrix
        F = randn(M,M - N);
    end
else
    if (realMat == 1)
        % Hadamard matrix;
        [A] = createCodingHadamardMat(MblockA, Nblock, numBlockLa, numBlockC, rpermA, Ja, noBlockErrorA);
        [F] = createCodingHadamardMat(MblockF, Nblock, numBlockL, numBlockC, rpermF, Jf, noBlockErrorF);
        if (showMat == 1); figure; imshow(F); figure; imshow(A'); pause; end
    end
end

disp([imageName, ' decoding']);
disp('true gamma'); disp(M / N);
disp('true alpha'); disp(sum(MblockF) / M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reconstructions

% lower the pixels value if required
if (strcmp(imageName, 'Lena') ); reduc = 256; else reduc = 1; end

% main loop
for ll = 1 : max(size(image) ) / sqrt(N)
    for cc = 1 : max(size(image) ) / sqrt(N)
        
        % patch of Lena
        S = double(image((ll - 1) * sqrt(N) + 1: ll * sqrt(N), (cc - 1) * sqrt(N) + 1 : cc * sqrt(N) ) ) ./ reduc;
        
        % Haar transform
        Sh = haar_transform(sqrt(N), levelHaar, S);
        
        % selection of the bigger Haar coefficients (in absolute value)
        [Ssort, IndSsort] = sort(abs(reshape(Sh, 1, N) ),'descend');
        SS = zeros(1, N); SS(IndSsort(1 : ceil(alphaWavelet * N) ) ) = Sh(IndSsort(1 : ceil(alphaWavelet * N) ) );
        
        % mixing of the components
        rp1 = randperm(N); S_mixed = SS(rp1);
        S_mixed = [S_mixed, zeros(1, sum(MblockA) - N) ];
        disp('number of concatenated zeros to the true signal'); disp(sum(MblockA) - N);
        
        % matrix null space
        if (hadamard == 0); A = null(F'); Y_perfect = A * S_mixed;
        else Y_perfect = MultSeededHadamardTransposeA(S_mixed, Ja, numBlockLa, numBlockC, MblockA, Nblock, rpermA, noBlockErrorA); end
        
        %                 A = randn(size(A));
        %                 Y_perfect = A' * S_mixed';
        
        % test of the operators
        if (realMat == 1)
            disp('***** beginning test of the operators, you should observe ~ 0 *****');
            max(Y_perfect - A' * S_mixed')
            SM = randn(M, 1);
            SN = randn(1, sum(MblockA) );
            test1 = MultSeededHadamardTransposeA(SN, Ja, numBlockLa, numBlockC, MblockA, Nblock, rpermA, noBlockErrorA);
            test_1 = SN * A;
            max(abs(test1' - test_1) )
            test2 = MultSeededHadamardTransposeSquarredA(SN, Ja, numBlockLa, numBlockC, MblockA, Nblock);
            test_2 = SN * A.^2;
            max(abs(test2' - test_2) )
            test3 = MultSeededHadamardA(SM, Ja, numBlockLa, numBlockC, MblockA, Nblock, rpermA, noBlockErrorA);
            test_3 = A * SM;
            max(abs(test3 - test_3) )
            test4 = MultSeededHadamardSquarredA(SM, Ja, numBlockLa, numBlockC, MblockA, Nblock);
            test_4 = A.^2 * SM;
            max(abs(test4 - test_4) )
        end
        
        % error due to noisy channel transmission
        k = floor(rhoCorrupt .* M);
        rp2 = randperm(M); error = [[sqrt(varBig) * randn(k, 1); zeros(M - k, 1) ] ]; error = error(rp2);
        small_error = randn(M, 1) * sqrt(varSmall);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Ideal decoding (no error due to transmission i.e perfect channel)
        
        if (recIdeal == 1)
            disp('ideal decoding')
            Y_id = Y_perfect + small_error;
            
            if (hadamard == 0);
                X_id(rp1) = A \ Y_id;
            else
                % algorithm properties
                My.signal_rho = sum(S_mixed ~= 0) / max(size(S_mixed) );
                My.signal = S_mixed;
                My.noBlockError = noBlockErrorA;
                My.M = N;
                My.N = M;
                My.J = Ja;
                My.numBlockL = numBlockLa;
                My.numBlockC = numBlockC;
                My.Mblock = MblockA;
                My.Nblock = Nblock;
                My.rp = rpermA;
                My.dump_mes = dumpMesA;
                My.nb_iter = nbIterIdeal;
                % operator or matrix based AMP method
                My.method = 'AMPhadamardSeededTransposeA';
                if (realMat == 1); My.method = 'AMP'; end
                % SparseGauss
                My.prior = 'SparseGauss';
                My.m_gauss = mean(S_mixed) / My.signal_rho;
                My.var_gauss = var(S_mixed) / My.signal_rho^2;
                if (realMat == 1); [resultsCS, n_and_e] = CSBP_Solver(Y_id, A', My);
                else [resultsCS, n_and_e] = CSBP_Solver(Y_id, [], My); end
                X_id(rp1) = resultsCS.av_mess(1 : N);
            end
            % inverse Haar transform
            X_id = haar_inverse(sqrt(N), levelHaar, reshape(X_id, sqrt(N), sqrt(N) ) );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Naive decoding (with no error correction)
        
        if (recNaive == 1)
            disp('naive decoding')
            My.nb_iter = nbIterNaive;
            if (hadamard == 0);
                X_naive(rp1) = A \ Y_err;
            else
                if (realMat == 1); [resultsCS, n_and_e] = CSBP_Solver(Y_err, A', My);
                else [resultsCS, n_and_e] = CSBP_Solver(Y_err, [], My); end
                X_naive(rp1) = resultsCS.av_mess(1 : N);
            end
            % inverse Haar transform
            X_naive = haar_inverse(sqrt(N), levelHaar, reshape(X_naive, sqrt(N), sqrt(N) ) );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% AMP decoding by reconstruction of the sparse error vector
        
        if (recAmp == 1)
            disp('AMP estimation of the error')
            % algorithm properties
            My.signal_rho = rhoCorrupt;
            My.signal = error + small_error;
            My.noBlockError = noBlockErrorF;
            My.N = M;
            My.M = sum(MblockF);
            My.J = Jf;
            My.numBlockL = numBlockL;
            My.numBlockC = numBlockC;
            My.Mblock = MblockF;
            My.Nblock = Nblock;
            My.rp = rpermF;
            My.dump_mes = dumpMesF;
            My.nb_iter = nbIterAmp;
            % operator or matrix based AMP method
            My.method = 'AMPhadamardSeeded';
            if (realMat == 1); My.method = 'AMP'; end
            % 2Gauss
            My.prior = '2Gauss';
            My.m_1_gauss = 0;
            My.m_2_gauss = 0;
            My.var_1_gauss = varSmall;
            My.var_2_gauss = varBig;
            if (realMat == 1); [results2G, n_and_e] = CSBP_Solver(Y_null', F', My);
            else [results2G, n_and_e] = CSBP_Solver(Y_null', [], My); end;
            Y_CS = Y_err - results2G.av_mess';
            
            disp('AMP decoding')
            if (hadamard == 0);
                X_CS(rp1) = A \ Y_CS;
            else
                % algorithm properties
                My.signal_rho = sum(S_mixed ~= 0) / max(size(S_mixed) );
                My.signal = S_mixed;
                My.noBlockError = noBlockErrorA;
                My.M = N;
                My.N = M;
                My.J = Ja;
                My.numBlockL = numBlockLa;
                My.numBlockC = numBlockC;
                My.Mblock = MblockA;
                My.Nblock = Nblock;
                My.rp = rpermA;
                My.dump_mes = dumpMesA;
                % operator or matrix based AMP method
                My.method = 'AMPhadamardSeededTransposeA';
                if (realMat == 1); My.method = 'AMP'; end
                % SparseGauss
                My.prior = 'SparseGauss';
                My.m_gauss = mean(S_mixed) / My.signal_rho;
                My.var_gauss = var(S_mixed) / My.signal_rho^2;
                if (realMat == 1); [resultsCS, n_and_e] = CSBP_Solver(Y_CS', A', My);
                else [resultsCS, n_and_e] = CSBP_Solver(Y_CS', [], My); end;
                X_CS(rp1) = resultsCS.av_mess(1 : N);
            end
            % inverse Haar transform
            X_CS = haar_inverse(sqrt(N), levelHaar, reshape(X_CS, sqrt(N), sqrt(N) ) );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% L1 decoding
        
        if (recL1 == 1)
            if (l1_AMP == 1)
                % AMP L1
                disp('L1 decoding')
                if (hadamard == 0);
                    X_L1(rp1) = A \ Y_CS;
                else
                    % algorithm properties
                    My.method = 'AMPhadamardSeededTransposeA';
                    if (realMat == 1); My.method = 'AMP'; end
                    My.prior = 'L1';
                    My.dump_mes = dumpMesA;
                    My.nb_iter = nbIterL1;
                    if (realMat == 1); [resultsCS, n_and_e] = CSBP_Solver(Y_CS', A', My);
                    else [resultsCS, n_and_e] = CSBP_Solver(Y_CS', [], My); end;
                    X_L1(rp1) = resultsCS.av_mess(1 : N);
                end
                % inverse Haar transform
                X_L1 = haar_inverse(sqrt(N), levelHaar, reshape(X_L1, sqrt(N), sqrt(N) ) );
            else
                % L1 optimization
                % Dantzig selector
                err_L1 = l1dantzig_pd(F * Y_null, F', [], Y_null, 3e-3, 1e-2, 1000);
                % regression step
                err_index = abs(err_L1) > sqrt(varSmall .* M); F_t = F';
                F_index_t = F_t(: , err_index);
                err_L1 = zeros(size(Y_err) ); err_L1(err_index) = F_index_t \ Y_null;
                Y_L1 = Y_err - err_L1;
                X_L1(rp1) = A \ Y_L1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        
        I_CS((ll - 1) * sqrt(N) + 1: ll * sqrt(N), (cc - 1) * sqrt(N) + 1: cc * sqrt(N) ) = reshape(X_CS, sqrt(N), sqrt(N) );
        I_L1((ll - 1) * sqrt(N) + 1: ll * sqrt(N), (cc - 1) * sqrt(N) + 1: cc * sqrt(N) ) = reshape(X_L1, sqrt(N), sqrt(N) );
        I_id((ll - 1) * sqrt(N) + 1: ll * sqrt(N), (cc - 1) * sqrt(N) + 1: cc * sqrt(N) ) = reshape(X_id, sqrt(N), sqrt(N) );
        I_naive((ll - 1) * sqrt(N) + 1: ll * sqrt(N), (cc - 1) * sqrt(N) + 1: cc * sqrt(N) ) = reshape(X_naive, sqrt(N), sqrt(N) );
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots

clf
subplot(3, 2, 1); imshow(image); xlabel('Original Image');
subplot(3, 2, 2); imshow(I_id); xlabel('Ideal');
subplot(3, 2, 3); imshow(I_CS); xlabel('AMP');
subplot(3, 2, 4); imshow(I_L1); xlabel('L1');
subplot(3, 2, 5); imshow(I_naive); xlabel('No error correction');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save

if (saveRes == 1); save(fileName,'I_id','I_CS','I_L1','I_naive','results_'); end