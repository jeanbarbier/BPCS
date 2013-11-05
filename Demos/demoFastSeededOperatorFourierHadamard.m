%% readme
% demo of fast, memory efficient and quasi-optimal CS reconstruction of a sparse-gauss real
% or complex signal by use of spatially-coupled structured operators.
% In the real signal case, each component is i.i.d with p(x) = (1 - rho) *
% delta(x) + rho * normal(x | mGauss, varGauss). In the complex case, the
% real and imaginary parts are i.i.d with same gaussian distribution normal(x | mGauss, varGauss).
% The final estimated signal is results.av_mess
% The bigger the number of blocks of the matrix, the slower the algorithm
% but a lower rate can be achieved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM AND PARAMETERS DEFINITION

clear

%% problem parameters
% type of the signal : 'real' (will use an Hadamard based operator) or 'complex' (will use a Fourier based operator)
type = 'complex';
% size of the signal (POWER OF 2 if type = 'real', will be improved soon...)
N = 5e5;
% ..and it's density/sparsity
rho = 0.3;
% mean and variance of the real signal or of the real and imaginary parts of the complex one
mGauss = 0;
varGauss = 1;
% global measurement rate (the true rate will be a little big bigger)
alphaGlobal = 0.5;
% coupling strenght
JJ = 0.2;
% number of blocks for the columns of the seeding matrix, it must divide N (POWER OF 2 if type = 'real')
numBlockC = 1;
% number of blocks for the rows of the seeding matrix
if (numBlockC == 1); numBlockL = 1; else numBlockL = numBlockC + 1; end
% measurement rate 1st block/seed, taken into acount if numBlockC > 1 : can be modified
if (numBlockC == 1); alphaCs(1) = alphaGlobal; else alphaCs(1) = 2 * rho; end
% measurement rate of the bulk blocks, taken into acount if numBlockC > 1 : not to be modified
if (numBlockC > 1); alphaCs(2 : numBlockL) = (numBlockL * alphaGlobal - alphaCs(1) ) / (numBlockL - 1); end
% coupling window (number of sub-diagonal blocks)
w = 6;

%% algorithm parameters
My = CSBP_Solver_Opt();
% max number of iterations
My.nb_iter = 1000;
% frequency of printing results
My.print = 1;
% convergence accuracy
My.conv = 1e-10;
% dumping
My.dump_mes = 0.;
% frequency of the plot of the block error between estimate and true signals (-1 for no use)
My.MSEbyBlock = 10; if (numBlockC == 1); My.MSEbyBlock = -1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO MODIFICATIONS AFTER THIS LINE ARE REQUIRED

%% creation of the signal
if (strcmp(type, 'complex') ); SnonZero = mGauss + (randn(1, floor(rho .* N) ) + 1i .* randn(1, floor(rho .* N) ) ) .* sqrt(varGauss);
else SnonZero = mGauss + randn(1, floor(rho .* N) ) .* sqrt(varGauss); end
S = zeros(1, N);
S(randperm(N, floor(rho .* N) ) ) = SnonZero;

%% matrix/operator block sizes
% number of variables per block in the seeded matrix, must divide N (POWER OF 2 if type = 'real')
Nblock = N / numBlockC;
% number of lines of the blocks in the seeded matrix
for (i = 1 : numBlockL) Mblock(i) = floor(alphaCs(i) * Nblock); end
disp('true measurement rate'); disp(sum(Mblock) / N);

%% block variance matrix
J = createSeededJ(numBlockL, numBlockC, JJ, w);

%% lines and signs randomization of the operators
[rp, flipedSigns] = createRandomLinesAndSignsPermutationForOperators(numBlockC, numBlockL, J, Mblock, Nblock);

%% algorithm inputs
if (strcmp(type, 'complex') ); My.prior = 'Complex'; My.method = 'AMPseededFourier'; else My.prior = 'SparseGauss'; My.method = 'AMPseededHadamard'; end
My.M = sum(Mblock);
My.N = N;
My.J = J;
My.numBlockL = numBlockL;
My.numBlockC = numBlockC;
My.Mblock = Mblock;
My.Nblock = Nblock;
My.rp = rp;
My.noBlockError = flipedSigns;
My.signal = S;
My.signal_rho = rho;
My.m_gauss = mGauss;
My.var_gauss = varGauss;
My.var_noise = 0;

%% measure
if (strcmp(type, 'complex') ); Y = MultSeededFourier(S, J, numBlockL, numBlockC, Mblock, Nblock, rp, flipedSigns);
else Y = MultSeededHadamard2(S, J, numBlockL, numBlockC, Mblock, Nblock, rp, flipedSigns); end

%% reconstruction
tic; [results, n_and_e, MSEt] = CSBP_Solver(Y, [], My); toc;