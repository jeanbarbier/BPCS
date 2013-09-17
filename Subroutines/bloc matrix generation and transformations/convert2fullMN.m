function [Jfull]=convert2fullMN(Jmat,Mvec,Nvec)
% Convert a small [L C] Jmat matrix to a full sparse [M N] one.

Jfull=sparse(sum(Mvec),sum(Nvec));
[L,C]=size(Jmat);
if (L~=max(size(Mvec))); error('error: L=length(Jmat) must be equal to length(Mvec)'); end;
if (C~=max(size(Nvec))); error('error: C=number_of_column(Jmat) must be equal to number_of_column(Mvec)'); end;

line=0;
for i=1:L
    line_prev=line+1; line=line+Mvec(i);
    column=0;
    for j=1:C
        column_prev=column+1; column=column+Nvec(j);
        Jfull(line_prev:line,column_prev:column)=Jmat(i,j);
    end
end

end

