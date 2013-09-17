clear all
% Parameters
Delta = 0; rho = 0.4; alpha = 0.5;
Emin = 1e-10; Emax = 1;
Vmin = 1e-10; Vmax = 1;
dE = 3;
dV = 3;
E = Emin; tt = 1;

while (E <= Emax)
    V = Vmin; t = 1;
    while (V <= Vmax)
        Frep(t,tt) = Frep_SparseBinary(rho,alpha,Delta,E,V);
        V_(1,t) = V;
        V = V .* dV;
        t = t + 1;
    end
    E_(1,tt) = E;
    E = E .* dE;
    tt = tt + 1;
end

Frep = Frep + abs(min(min(Frep(:,:) ) ) );
surf(E_,V_,Frep); figure(gcf);
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlabel('E')
ylabel('V')

figure;
% Plot on the Nishimori line (only if Frep is a squarre matrix)
loglog(E_,diag(Frep) );
xlabel('E = V')
ylabel('Bethe F')