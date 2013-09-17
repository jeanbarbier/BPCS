function [S] = S_SparseBinary(sze,rho)
% Creates a binary (0,1) signal of density rho

S = (rand(1,sze) < rho);

end

