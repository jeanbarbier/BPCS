clear

% generation of the path
addpath(genpath('../') );
% proportion of measurements
alpha = 0.2;
% size of the signal..
N = 1e6; M = ceil(alpha .* N);
% ..and it's density
rho = 0.1;
% creation of the signal
mGauss = 0;
varGauss = 1;
SnonZero = mGauss + (randn(1, floor(rho .* N) ) + i .* randn(1, floor(rho .* N) ) ) .* sqrt(varGauss);
S = zeros(1, N);
S(randperm(N, floor(rho .* N) ) ) = SnonZero;
% coupling strenght
JJ = 0.2;
% coupling window
w = 2;

% properties of the measurement matrix
numBlockC = 10;
numBlockL = 10;
Nblock = N ./ numBlockC;
J = createSeededHadamardJ(numBlockL, numBlockC, JJ, w);
for l = 1 : numBlockL
    Mblock(l) = M ./ numBlockL;
    %     noBlockError{l} = randperm(Mblock(l), floor(Mblock(l) ./ 2) );
    noBlockError{l} = [];
    for c = 1 : numBlockC; rp{l, c} = randperm(Nblock - 1) + 1; end
end

% measure
Y = MultSeededFourier(S, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError);
% fftmtx = dftmtx(N);
% G = fftmtx(rp{1,1}, :);
% Y = G * S';

% algorithm properties
My = CSBP_Solver_Opt();
My.nb_iter = 100;
My.save_speed = 0;
My.save_memory = 1 - My.save_speed;
My.print = 1;
My.conv = 1e-10;
My.learn = 0;
My.signal_rho = rho;
My.var_noise = 0;
My.option_noise = 0;
My.remove_mean = 0;
My.signal = S;
My.dump_mes = 0.5;
My.prior = 'Complex';
My.method = 'AMPseededFourier';
% My.method = 'AMPcomplex';
My.M = M;
My.N = N;
My.J = J;
My.numBlockL = numBlockL;
My.numBlockC = numBlockC;
My.Mblock = Mblock;
My.Nblock = Nblock;
My.rp = rp;
My.noBlockError = noBlockError;
My.m_gauss = mGauss;
My.var_gauss = varGauss;

% reconstruction
if (strcmp(My.method,'AMPseededFourier') ); tic; [results,n_and_e] = CSBP_Solver(Y, [], My); toc;
else tic; [results,n_and_e] = CSBP_Solver(Y, G, My); toc; end