function [E_new_rep,V_new_rep,E_new_cav,V_new_cav] = DensEvoSparseGauss(rep,cav,rho,alpha,Delta,E,V)
% Density evolution iteration for a gaussian signal with zero mean and unit variance.
minn = -10; maxx = -minn;
E_new_rep = 0; E_new_cav = 0; V_new_rep = 0; V_new_cav = 0;

if (cav == 1)
    Z = @(S2,R,rho) (1 - rho) ./ sqrt(2 .* pi .* S2) .* exp(-R.^2 ./ (2 .* S2) ) + rho ./ sqrt(2 .* pi .* (S2 + 1) ) .* exp(-0.5 .* R.^2 ./ (1 + S2) );
    f_a = @(S2,R,rho) (rho ./ sqrt(2 .* pi .* (S2 + 1) ) .* exp(-0.5 .* R.^2 ./ (1 + S2) ) .* R ./ (S2 + 1) ) ./ Z(S2,R,rho);
    f_b = @(S2,R,rho) (rho ./ sqrt(2 .* pi .* S2) .* (S2.^(-1) + 1).^(-3 ./ 2) .* exp(-0.5 .* R.^2 ./ (1 + S2) ) .* (1 + (R ./ S2).^2 ./ (1 ./ S2 + 1) ) ) ./ Z(S2,R,rho);
    f_c = @(S2,R,rho) max(1e-50, f_b(S2,R,rho) - f_a(S2,R,rho).^2);
    G = @(z) exp(-0.5 .* z.^2) ./ sqrt(2 .* pi);
    
    int_ = integral(@(s) G(s) .* f_c((Delta + V) ./ alpha,s .* sqrt((E + Delta) ./ alpha + 1),rho),minn,maxx,'AbsTol',1e-20);
    V_new_cav = (1 - rho) .* integral(@(s) G(s) .* f_c((Delta + V) ./ alpha,s .* sqrt((E + Delta) ./ alpha),rho),minn,maxx,'AbsTol',1e-20) + rho .* int_;
    E_new_cav = (1 - rho) .* integral(@(s) G(s) .* f_a((Delta + V) ./ alpha,s .* sqrt((E + Delta) ./ alpha),rho).^2,minn,maxx,'AbsTol',1e-20) + rho .* (1 - 2 .* alpha ./ (Delta + V) .* int_ + integral(@(s) G(s) .* f_a((Delta + V) ./ alpha,s .* sqrt((E + Delta) ./ alpha + 1),rho).^2,minn,maxx,'AbsTol',1e-20) );
end

if (rep == 1)
    G = @(z) exp(-0.5 .* z.^2) ./ sqrt(2 .* pi);
    h1 = @(rho,alpha,Delta,E,V,z) (1 - rho) .* exp(-z.^2 .* alpha .* (Delta + E) ./ (2.* (Delta + V) .* (Delta + V + alpha) ) ) + rho .* sqrt((Delta + V) ./ (Delta + V + alpha) );
    h2 = @(rho,alpha,Delta,E,V,z) (1 - rho) .* exp(-z.^2 .* alpha .* (Delta + E + alpha) ./ (2.* (Delta + V) .* (Delta + V + alpha) ) ) + rho .* sqrt((Delta + V) ./ (Delta + V + alpha) );
    dE_h1 = @(rho,alpha,Delta,E,V,z) (rho - 1) .*  z.^2 .* alpha ./ (2.* (Delta + V) .* (Delta + V + alpha) ) .* exp(-z.^2 .* alpha .* (Delta + E) ./ (2.* (Delta + V) .* (Delta + V + alpha) ) );
    dV_h1 = @(rho,alpha,Delta,E,V,z) (1 - rho) .*  z.^2 .* alpha .* (Delta + E) ./ ((Delta + V) .* (Delta + V + alpha) ).^2 .* exp(-z.^2 .* alpha .* (Delta + E) ./ (2.* (Delta + V) .* (Delta + V + alpha) ) ) .* (Delta + V + 0.5 .* alpha) + 0.5 .* rho .* sqrt(1 + alpha ./ (Delta + V) ) .* alpha ./ (Delta + V + alpha).^2;
    dE_h2 = @(rho,alpha,Delta,E,V,z) (rho - 1) .*  z.^2 .* alpha ./ (2.* (Delta + V) .* (Delta + V + alpha) ) .* exp(-z.^2 .* alpha .* (Delta + E + alpha) ./ (2.* (Delta + V) .* (Delta + V + alpha) ) );
    dV_h2 = @(rho,alpha,Delta,E,V,z) (1 - rho) .*  z.^2 .* alpha .* (Delta + E + alpha) ./ ((Delta + V) .* (Delta + V + alpha) ).^2 .* exp(-z.^2 .* alpha .* (Delta + E + alpha) ./ (2.* (Delta + V) .* (Delta + V + alpha) ) ) .* (Delta + V + 0.5 .* alpha) + 0.5 .* rho .* sqrt(1 + alpha ./ (Delta + V) ) .* alpha ./ (Delta + V + alpha).^2;
    dE_f = -alpha ./ 2 .* V ./ (Delta + V).^2;
    dV_f = -alpha ./ (2 .* (Delta + V).^3) .* ((V - E) .* (V - Delta) - rho .* (Delta + V) );
    dE_g = alpha ./ (2 .* (Delta + V) .* (Delta + V + alpha) );
    dV_g = -alpha .* (Delta + E + rho .* alpha) .* (Delta + V + 0.5 .* alpha) ./ ((Delta + V) .* (Delta + V + alpha) ).^2;
    
    dE_Frep = dE_f + dE_g + (1 - rho) .* integral(@(s) G(s) .* dE_h1(rho,alpha,Delta,E,V,s) ./ h1(rho,alpha,Delta,E,V,s),minn,maxx,'AbsTol',1e-20) + rho .* integral(@(s) G(s) .* dE_h2(rho,alpha,Delta,E,V,s) ./ h2(rho,alpha,Delta,E,V,s),minn,maxx,'AbsTol',1e-20);
    dV_Frep = dV_f + dV_g + (1 - rho) .* integral(@(s) G(s) .* dV_h1(rho,alpha,Delta,E,V,s) ./ h1(rho,alpha,Delta,E,V,s),minn,maxx,'AbsTol',1e-20) + rho .* integral(@(s) G(s) .* dV_h2(rho,alpha,Delta,E,V,s) ./ h2(rho,alpha,Delta,E,V,s),minn,maxx,'AbsTol',1e-20);
    V_new_rep = V + dV_Frep;
    E_new_rep = E + dE_Frep;
end


end

