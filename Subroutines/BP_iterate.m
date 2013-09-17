%disp('BP_iterate');
%1.) This function always remove mean.
%2.) It always assume that Gmean is not given, which should be computed by
%the algorithm, so dimension of Gmean is [M, N]!
assert(Gmean_dimention>1);

%compute vmu omega 
%version 1, matrix multiplication, much more useless computations are paid.
%vmu = diag(Geff2*F_b_c); % diag([M,N]*[N,M])  = [M,1]
%omega = diag(Geff*F_a_c); %[M,1]
%version 2, loop, should be also not efficient in Matlab.
%for i=1:M
%    vmu(i)=Geff2(i,:)*F_b_c(:,i);%[1,N]*[N,1]
%    omega(i)=Geff(i,:)*F_a_c(:,i);%[1,N]*[N,1]
%end
%version 3, fastest one
%vmu = sum(Geff2'.*F_b_c)';
%omega=sum(Geff'.*F_a_c)';
%version 4, use sparsity
vmu = (sum(G2'.*F_b_c-2*GGm'.*F_b_c)+Gm1d2*F_b_c)';
omega=(sum(G'.*F_a_c)-Gmean1d*F_a_c)';
%version 1
%vmu2i=vmu*ones(1,N)-Geff2.*F_b_c'; %[M,N]
%omegamu2i=omega*ones(1,N)-Geff.*F_a_c'; %[M,N]'
%version 2
vmu2i=bsxfun(@minus,vmu,Geff2.*F_b_c');
omegamu2i=bsxfun(@minus,omega,Geff.*F_a_c');

%compute Amu2i and Bmu2i
Amu2i=Geff2./(var_noise+vmu2i);
Bmu2i=Geff.*(Yeffext-omegamu2i)./(var_noise+vmu2i);

%compute C and D
%version 1
%C=Amu2i'*ones(M,1);%[N,1]
%D=Bmu2i'*ones(M,1);%[N,1]
%Ci2mu=(C*ones(1,M))-Amu2i'; %?
%Di2mu=(D*ones(1,M))-Bmu2i'; %?
%version 2
C=sum(Amu2i)';%[N,1]
D=sum(Bmu2i)';%[N,1]
Ci2mu=bsxfun(@minus,C,Amu2i');
Di2mu=bsxfun(@minus,D,Bmu2i');

%compute F_a_c, F_b_c, F_a and F_b
S2_c_new=1./Ci2mu;%[N,M]
R_c_new=Di2mu./Ci2mu;%[N,M]
S2_new=1./C'; %[1,N]
R_new=D'./C'; %[1,N]


%damping 
if (t==1)
    R=R_new;
    S2=S2_new;
    S2_c=S2_c_new;
    R_c=R_c_new;
else
    R=damp_mes*R_new+(1-damp_mes)*R;
    S2=damp_mes*S2_new+(1-damp_mes)*S2;
    R_c=damp_mes*R_c_new+(1-damp_mes)*R_c;
    S2_c=damp_mes*S2_c_new+(1-damp_mes)*S2_c;
end

F_a_c=Fun_a(S2_c,R_c,rho,exp_gauss,sg2);
F_b_c=Fun_b(S2_c,R_c,rho,exp_gauss,sg2);
F_a_new=Fun_a(S2,R,rho,exp_gauss,sg2);
F_b=Fun_b(S2,R,rho,exp_gauss,sg2);
