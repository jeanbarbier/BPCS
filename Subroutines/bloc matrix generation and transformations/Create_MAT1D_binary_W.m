function [F,Jmat,Mvec,Nvec] = Create_MAT1D_binary_W(N,L,J,alpha1,alpha2,W)

N2=floor(N/L);
N1=N-(L-1)*N2;
M1=floor(N1*alpha1);
M2=floor(N2*alpha2);
M=M1+(L-1)*M2;
alpha=M/N

[M N M1 N1 M2 N2]

%Matrix J
Jmat=sparse(L,L);
for w=0:W
    for i=1:L
        j=i-w;
        if j>0
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
F1=2*(rand(M1,N1)>0.5)-1;
F2=(2*(rand(M1,N2)>0.5)-1).*(rand(M1,N2)<J);
F=[F1 F2 sparse(M1,N-N1-N2)];
if(L>1)
%All lines (expect the first and last)
for i=2:L
 F1=sparse(0,0);

 if Jmat(i,1)==0
      F1=[F1 sparse(M2,N1)];
 else
      F1=[F1 (2*(rand(M2,N1)>0.5)-1).*(rand(M2,N1)<Jmat(i,1))];
 end
 
 for j=2:L
    if Jmat(i,j)==0
      F1=[F1 sparse(M2,N2)];
    else
      F1=[F1 (2*(rand(M2,N2)>0.5)-1).*(rand(M2,N2)<Jmat(i,j))];
    end
 end
 F=[F; F1];
end

end
