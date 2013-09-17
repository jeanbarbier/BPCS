function [S] = S_SparseConstant(size,rho,down,up)
% Creates a sparse bounded constant signal of density rho in [down,up].

SS = [zeros(1,size - ceil(rho .* size)),(up - down) .* rand(1,ceil(rho .* size) ) + down];
S = (intrlv(SS,randperm(size) ) );

end

