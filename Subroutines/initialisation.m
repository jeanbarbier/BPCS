if (strcmpi(opt.method,'AMP') || (strcmpi(opt.method,'AMPhb') ) || (strcmpi(opt.method,'AMPcomplex') ) || (strcmpi(opt.method,'AMPhadamard') ) || (strcmpi(opt.method,'AMPseededHadamard') ) || (strcmpi(opt.method,'AMPseededHadamardTranspose') ) || (strcmpi(opt.method,'AMPseededFourier') ) )
    if (opt.remove_mean == 0); W = Y; else W = Yeff; end; V = ones(1,M);
    R_init = zeros(1,N); S2_init = zeros(1,N); av_mess_init = zeros(1,N); var_mess_init = opt.signal_rho .* ones(1,N);
end

if (strcmpi(opt.method,'AMPseededHadamardTransposeA') )
    if (opt.remove_mean == 0); W = Y; else W = Yeff; end; V = ones(1,N);
    R_init = zeros(1,M); S2_init = zeros(1,M); av_mess_init = zeros(1,M); var_mess_init = opt.signal_rho .* ones(1,M);
end

if (strcmpi(opt.method,'AMPh') )
    if (opt.remove_mean == 0); W = Y; else W = Yeff; end; V = 1;
    R_init = zeros(1,N); S2_init = zeros(1,N); av_mess_init = zeros(1,N); var_mess_init = opt.signal_rho .* ones(1,N);
end

if (strcmp(opt.method,'BP') )
    small_value=0.0001;
    [idxi,idxj]=find(G);
    F_a_c=sparse(idxi,idxj,small_value,M,N); F_a_c=F_a_c.'; %cavity average, [N,M] with sparse representation.
    F_b_c=sparse(idxi,idxj,rho,M,N); F_b_c=F_b_c.'; %cavity variance, [N,M] with sparse representation.
    F_a=zeros(1,N);
    F_b=ones(1,N)*rho;
    C=zeros(N,1); %marginal message for i [N,1]
    D=zeros(N,1); %marginal message for i [N,1]
    Amu2i=sparse(idxi,idxj,small_value,M,N); %cavity message from mu to i, [M,N]
    Bmu2i=sparse(idxi,idxj,small_value,M,N); %cavity message from mu to i, [M,N]
    vmu=zeros(M,1); %marginal message for mu [M,1]
    omega=zeros(M,1); %marginal message for mu [M,1]
    if (opt.remove_mean == 0); W = Y; else W = Yeff; end; V = ones(1,M);
end

free_entropy = [];
