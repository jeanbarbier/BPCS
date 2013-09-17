% Initialisation 
%disp('BP_init');
small_value=0.0001;

[idxi,idxj]=find(G);
F_a_c=sparse(idxi,idxj,small_value,M,N); F_a_c=F_a_c'; %cavity average, [N,M] with sparse representation.
F_b_c=sparse(idxi,idxj,rho,M,N); F_b_c=F_b_c'; %cavity variance, [N,M] with sparse representation.
F_a=zeros(1,N);
F_b=ones(1,N)*rho;
C=zeros(N,1); %marginal message for i [N,1]
D=zeros(N,1); %marginal message for i [N,1]

Amu2i=sparse(idxi,idxj,small_value,M,N); %cavity message from mu to i, [M,N]
Bmu2i=sparse(idxi,idxj,small_value,M,N); %cavity message from mu to i, [M,N]
vmu=zeros(M,1); %marginal message for mu [M,1]
omega=zeros(M,1); %marginal message for mu [M,1]

if (remove_mean==0)
    W=Y;
else
    W=Y_eff;
end

disp('Done');