function [F] = Frep_SparseGauss(rho,alpha,Delta,E,V)
% Replica free energy for a gaussian signal with zero mean and unit variance.

f = rho .* sqrt((V + Delta) ./ (V + alpha + Delta) );
int_1 = @(z) exp(-z.^2 ./ 2) .* log((1 - rho) .* exp(-0.5 .* z.^2 .* alpha .* (Delta + E) ./ ((Delta + V) .* (Delta + V + alpha) ) ) + f);
int_2 = @(z) exp(-z.^2 ./ 2) .* log((1 - rho) .* exp(-0.5 .* z.^2 .* alpha .* (Delta + E + alpha) ./ ((Delta + V) .* (Delta + V + alpha) ) ) + f);
int1 = integral(int_1,-10,10); 
int2 = integral(int_2,-10,10);
part1 = -0.5 .* alpha .* (Delta ./ (Delta + V) + log(Delta + V) + rho ./ (Delta + V) - V .* (V - E) ./ (Delta + V).^2 );
part2 = 0.5 .* alpha ./ ((Delta + V) .* (Delta + V + alpha) ) .* ((1 - rho) .* (Delta + E) + rho .* (Delta + E + alpha) );
part3 = (1 - rho) ./ sqrt(2 .* pi) .* int1;
part4 = rho ./ sqrt(2 .* pi) .* int2;
F = part1 + part2 + part3 + part4; 

end

