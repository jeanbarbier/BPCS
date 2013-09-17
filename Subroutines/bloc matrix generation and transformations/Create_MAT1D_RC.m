function [F,Jmat,Mvec,Nvec] = Create_MAT1D_RC(N,M,L,J1,J2,alpha2)
% SYNTAX:
% [F,Jmat,Mvec,Nvec] = Create_MAT1D(N,L,J1,J2,alpha1,alpha2)

N2=floor(N/L);
N1=N-(L-1)*N2;
M2=floor(N2*alpha2);
M1 = M - (L - 1) .* M2;
disp(M/N)
[M N M1 N1 M2 N2 alpha2 M1 ./ N1]

%Matrix J
Jmat=sparse(L,L);
for i=1:L
    Jmat(i,i)=1;
    if (i<L)
        Jmat(i,i+1)=J2;
    end
    if (i<L)
        Jmat(i+1,i)=J1;
    end
end

%Matrix M
Mvec=zeros(1,L);
Mvec(1)=M1;
for i=2:L
    Mvec(i)=M2;
end


%Matrix N
Nvec=zeros(1,L);
Nvec(1)=N1;
for i=2:L
    Nvec(i)=N2;
end

%Matrix F (sparse!)
%First line
F1=randn(M1,N1);
F2=randn(M1,N2)*sqrt(J2);
F=[F1 F2 sparse(M1,N-N1-N2)];
%Second line
F1=randn(M2,N1)*sqrt(J1);
F2=randn(M2,N2);
if (L>2)
    F3=randn(M2,N2)*sqrt(J2);
    F=[F; F1 F2 F3 sparse(M2,N-N1-2*N2)];
else
    F=[F; F1 F2];
end


if(L>2)
    %All lines (expect the first second and last)
    for i=3:L-1
        F1=randn(M2,N2)*sqrt(J1);
        F2=randn(M2,N2);
        F3=randn(M2,N2)*sqrt(J2);
        F=[F; sparse(M2,N1+(i-3)*N2) F1 F2 F3 sparse(M2,N-N1-i*N2)];
    end
    
    %last line
    F1=randn(M2,N2)*sqrt(J1);
    F2=randn(M2,N2);
    F=[F; sparse(M2,N1+(L-3)*N2) F1 F2 ];
end

end


