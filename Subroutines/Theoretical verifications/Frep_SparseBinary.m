function [F] = Frep_SparseBinary(rho,alpha,Delta,E,V)
% Replica free energy for a binary signal.

int_1 = @(z) exp(-z.^2 ./ 2) .* log(rho .* exp(-0.5 .* alpha ./ (Delta + V) + z .* sqrt(alpha .* (Delta + E) ) ./ (Delta + V) ) + (1 - rho) );
int_2 = @(z) exp(-z.^2 ./ 2) .* log(rho + (1 - rho) .* exp(-0.5 .* alpha ./ (Delta + V) - z .* sqrt(alpha .* (Delta + E) ) ./ (Delta + V) ) );
int1 = quadgk(int_1,-10,10); 
int2 = quadgk(int_2,-10,10);
part1 = -0.5 .* alpha .* ((Delta + E) ./ (Delta + V) + log(Delta + V) + 0.5 .* (rho.^2 - E) ./ (Delta + V) - V ./ (Delta + V) + V .* (E + Delta) ./ (Delta + V).^2);
part2 = (1 - rho) ./ sqrt(2 .* pi) .* int1;
part3 = rho ./ sqrt(2 .* pi) .* int2 + rho .* 0.5 .* alpha ./ (Delta + V);
F = part1 + part2 + part3; 

end

