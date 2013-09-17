function [Fmean,Ymean] = Remove_Mean_GT(Mvec,F,Y)

L=max(size(Mvec));

%Matrix Fmean (sparse!)
Fmean=[];
Ymean=[];
Mstart=1;Mend=Mvec(1);
for i=1:L
    Msize=1+Mend-Mstart;
    
    ThisFmean=sparse(sum(F(Mstart:Mend,:)))/(Msize);
    ThisYmean=sum(Y(Mstart:Mend))/(Msize);
    

    Fmean=[Fmean;ones(Msize,1)*ThisFmean];
    Ymean=[Ymean;ones(Msize,1)*ThisYmean];
    

    Mstart=Mend+1;
    if (i<L)
        Mend=Mend+Mvec(i+1);
    end
end




