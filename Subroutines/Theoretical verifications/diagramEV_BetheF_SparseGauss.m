clear all
% Parameters
rep = 1; cav = 1;
Delta = 0; rho = 0.33; alpha = 0.5;
Emin = 1e-2; Emax = 0.12;
Vmin = 1e-2; Vmax = 0.12;
dE = 0.003; dV = dE;
E = Emin; tt = 1;
disp((Emax - Emin) .* (Vmax - Vmin) ./ (dE .* dV) )

while (E <= Emax)
    
    V = Vmin; t = 1;
    while (V <= Vmax)
        % Free entropy
        Frep_ = Frep_SparseGauss(rho,alpha,Delta,E,V); Frep(tt,t) = Frep_;
        
        % Density evolution
        [E_new_rep,V_new_rep,E_new_cav,V_new_cav] = DensEvoSparseGauss(rep,cav,rho,alpha,Delta,E,V);
        if (cav == 1); DE_E_cav(tt,t) = (E_new_cav - E) ./ sqrt((V - V_new_cav).^2 + (E - E_new_cav).^2); DE_V_cav(tt,t) = (V_new_cav - V) ./ sqrt((V - V_new_cav).^2 + (E - E_new_cav).^2); end;
        if (rep == 1); DE_E_rep(tt,t) = (E_new_rep - E) ./ sqrt((V - V_new_rep).^2 + (E - E_new_rep).^2); DE_V_rep(tt,t) = (V_new_rep - V) ./ sqrt((V - V_new_rep).^2 + (E - E_new_rep).^2); end;
        
        % New value of Frep to see if it max or min
        if (cav == 1); Frep_new_cav = Frep_SparseGauss(rho,alpha,Delta,E_new_cav,V_new_cav); end;
        if (rep == 1); Frep_new_rep = Frep_SparseGauss(rho,alpha,Delta,E_new_rep,V_new_rep); end;
        if (cav == 1); if (Frep_new_cav > Frep_); DIR_cav(tt,t) = 1; else DIR_cav(tt,t) = -1; end; end;
        if (rep == 1); if (Frep_new_rep > Frep_); DIR_rep(tt,t) = 1; else DIR_rep(tt,t) = -1; end; end;
        
        V_(tt,t) = V; E_(tt,t) = E; V = V + dV; t = t + 1;
    end
    
    E = E + dE;  tt = tt + 1;
end

% % 3D plot of F
% % Frep = Frep + abs(min(min(Frep(:,:) ) ) );
% figure('Name','Frep');
% surf(E_(:,1),V_(1,:),Frep');
% xlabel('E'); ylabel('V');
% % set(gca, 'xscale', 'log');
% % set(gca, 'yscale', 'log');
% % set(gca, 'zscale', 'log');

% figure;
% Fflo = Fflo + abs(min(min(Fflo(:,:) ) ) );
% figure('Name','Fflo');
% surf(E_(:,1),V_(1,:),Fflo');
% xlabel('E'); ylabel('V');
% set(gca, 'xscale', 'log');
% set(gca, 'yscale', 'log');
% set(gca, 'zscale', 'log');

% 3D plot of the derivative wrt E and V
for (i = 1 : tt - 2)
    for (j = 1 : t - 2)
        dE_Frep(i,j) = (Frep(i + 1,j) - Frep(i,j) ) ./ dE;
        dV_Frep(i,j) = (Frep(i,j + 1) - Frep(i,j) ) ./ dV;
    end
end

% cavity density evolution flow + derivative of Frep
if (cav == 1)
figure('Name','dE_Frep + DE flow_cav');
% contourf(E_(1 : end - 1,1),V_(1,1 : end - 1),(dE_Frep' < 0) );
hold on;
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == 1),DE_V_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == 1),'g');
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == -1),DE_V_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == -1),'r');
% hold off;
xlabel('E'); ylabel('V');

figure('Name','dV_Frep + DE flow_cav');
% contourf(E_(1 : end - 1,1),V_(1,1 : end - 1),(dV_Frep' < 0) );
hold on;
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == 1),DE_V_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == 1),'g');
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == -1),DE_V_cav(1 : end - 1,1 : end - 1) .* (DIR_cav(1 : end - 1,1 : end - 1) == -1),'r');
% hold off;
xlabel('E'); ylabel('V');
end

% replica density evolution flow + derivative of Frep
if (rep == 1)
figure('Name','dE_Frep + DE flow_rep');
% contourf(E_(1 : end - 1,1),V_(1,1 : end - 1),(dE_Frep' < 0) );
hold on;
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == 1),DE_V_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == 1),'g');
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == -1),DE_V_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == -1),'r');
hold off;
xlabel('E'); ylabel('V');

figure('Name','dV_Frep + DE flow_rep');
% contourf(E_(1 : end - 1,1),V_(1,1 : end - 1),(dV_Frep' < 0) );
hold on;
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == 1),DE_V_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == 1),'g');
quiver(E_(1 : end - 1,1 : end - 1),V_(1 : end - 1,1 : end - 1),DE_E_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == -1),DE_V_rep(1 : end - 1,1 : end - 1) .* (DIR_rep(1 : end - 1,1 : end - 1) == -1),'r');
hold off;
xlabel('E'); ylabel('V');
end;

% % Plot on the Nishimori line (only if Frep is a squarre matrix)
% figure;
% semilogx(E_,diag(Frep) );
% xlabel('E = V'); ylabel('Bethe F');

% % Plot on the direction perpendicular to the Nishimori line (only if Frep is a squarre matrix)
% figure;
% semilogx(E_,diag(fliplr(Frep) ) );
% xlabel('perpendicular to Nishimori')
% ylabel('Bethe F')

% % Plot the flow of the algorithm given by density evolution with information on the increase or decrease of F
% figure('Name','Density Evolution Flow');
% quiver(E_,V_,DE_E .* (DIR == 1),DE_V .* (DIR == 1),'g');
% hold on
% quiver(E_,V_,DE_E .* (DIR == -1),DE_V .* (DIR == -1),'r');
% xlabel('E'); ylabel('V');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONLY NISHIMORI:

% clear all
% % Parameters
% Delta = 0; rho = 0.35; alpha = 0.48;
% Emin = 1e-5; Emax = 1;
% dE = 1.1;
% E = Emin; t = 1;
%
% while (E <= Emax)
%     V = E;
%     Frep(1,t) = Frep_SparseGauss(rho,alpha,Delta,E,V);
%     E_(1,t) = E;
%     E = E .* dE;
%     t = t + 1;
% end
% semilogx(E_,Frep);
% xlabel('E')