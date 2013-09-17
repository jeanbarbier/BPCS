clear all
% Computation parameters and initialisation
epsilon = 1e-14; dalpha = 5e-5; min = -20; max = -min; dump = 0.9; alpha = 1e-5; t = 1; alpha_max = 0.01;
% Signal parameters (here, for a Laplace prior)
beta = 1; Delta = 0; s2m = 2 ./ beta.^2;

% Laplace prior
Ps = @(s,beta) beta ./ 2 .* exp(-beta .* abs(s) );

uu = @(s,alpha,V,E,Delta) (Delta + V + s .* alpha) ./ sqrt(alpha .* (Delta + E) );
vv = @(s,alpha,V,E,Delta) (Delta + V - s .* alpha) ./ sqrt(alpha .* (Delta + E) );
integrandV = @(s,alpha,V,E,Delta,beta) Ps(s,beta) .* (erfc(uu(s,alpha,V,E,Delta) ./ sqrt(2) ) + erfc(vv(s,alpha,V,E,Delta) ./ sqrt(2) ) );
integrandE1 = @(s,alpha,V,E,Delta,beta) Ps(s,beta) .* (erfc(uu(s,alpha,V,E,Delta) ./ sqrt(2) ) + erfc(vv(s,alpha,V,E,Delta) ./ sqrt(2) ) );
integrandE2 = @(s,alpha,V,E,Delta,beta) Ps(s,beta) .* (exp(-uu(s,alpha,V,E,Delta).^2 ./ 2 ) + exp(-vv(s,alpha,V,E,Delta).^2 ./ 2 ) );
integrandE3 = @(s,alpha,V,E,Delta,beta) Ps(s,beta) .* (exp(-uu(s,alpha,V,E,Delta).^2 ./ 2 ) .* uu(s,alpha,V,E,Delta) ./ sqrt(2 .* pi) + exp(-vv(s,alpha,V,E,Delta).^2 ./ 2 ) .* vv(s,alpha,V,E,Delta) ./ sqrt(2 .* pi) + 0.5 .* (erfc(uu(s,alpha,V,E,Delta) ./ sqrt(2) ) + erfc(vv(s,alpha,V,E,Delta) ./ sqrt(2) ) ) );
integrandE4 = @(s,alpha,V,E,Delta,beta) Ps(s,beta) .* s.^2 .* (erfc(uu(s,alpha,V,E,Delta) ./ sqrt(2) ) + erfc(vv(s,alpha,V,E,Delta) ./ sqrt(2) ) );
integrandE = @(s,alpha,V,E,Delta,beta) ((Delta + V) ./ alpha).^2 ./ 2 .* integrandE1(s,alpha,V,E,Delta,beta) - (Delta + V) ./ alpha .* sqrt(2 .* (E + Delta) ./ (alpha .* pi) ) .* integrandE2(s,alpha,V,E,Delta,beta) + (Delta + E) ./ alpha  .* integrandE3(s,alpha,V,E,Delta,beta) - 0.5 .* integrandE4(s,alpha,V,E,Delta,beta);

while (alpha <= alpha_max)
    
    V = 100; E = V; dE = 100; dV = dE; alpha
    
    while ((dV > epsilon) || (dE > epsilon) )
        V_new = 0.5 .* (Delta + V) ./ alpha .* integral(@(s)integrandV(s,alpha,V,E,Delta,beta),min,max);
        E_new = s2m + integral(@(s)integrandE(s,alpha,V,E,Delta,beta),min,max);
        
        dE = abs(E - E_new); dV = abs(V - V_new);
        
        V = dumping(V,V_new,dump); E = E_new;
    end
    
    alpha_(1,t) = alpha; E_(1,t) = E;  V_(1,t) = V;
    alpha = alpha + dalpha; t = t + 1;
    
end



% V = 100; E = V; dE = 100; dV = dE; alpha = 0.9; dump = 0.9;
%
% while ((dV > epsilon) || (dE > epsilon) )
%     V_new = 0.5 .* (Delta + V) ./ alpha .* integral(@(s)integrandV(s,alpha,V,E,Delta,beta),min,max);
%     E_new = s2m + integral(@(s)integrandE(s,alpha,V,E,Delta,beta),min,max);
%
%     dE = abs(E - E_new); dV = abs(V - V_new);
%
%     V = dumping(V,V_new,dump); E = E_new;
% end
%
% [E,V]