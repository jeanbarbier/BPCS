function [F,Jmat,Mvec,Nvec] = Create_MAT1D_gauss(N,L,J,alpha1,alpha2)

N2=floor(N/L);
N1=N-(L-1)*N2;
M1=floor(N1*alpha1);
M2=floor(N2*alpha2);
M=M1+(L-1)*M2;
alpha=M/N

[M N M1 N1 M2 N2]

%Matrix J
Jmat=sparse(L,L);
for i=1:L
    for j=1:L
        if (i>=j)
          Jmat(i,j)=1;
        end
    end
end
    
for i=1:L
  if (i<L)
    Jmat(i,i+1)=J;
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
F2=randn(M1,N2)*sqrt(J);
F=[F1 F2 sparse(M1,N-N1-N2)];
if(L>1)
%All lines (expect the first and last)
for i=2:L-1
  F1=randn(M2,N2*(i-1)+N1);
  F2=randn(M2,N2)*sqrt(J);
  F=[F; F1 F2 sparse(M2,N-N1-i*N2)];
end

%last line
F1=randn(M2,N);
F=[F; F1];
end
