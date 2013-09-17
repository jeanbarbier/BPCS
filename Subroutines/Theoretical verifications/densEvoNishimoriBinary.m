% plot of MSE(alpha) via density evolution for a binary signal [1, -1] distributed as p(x) = rho *
% delta(x - 1) + (1 - rho) * delta(x + 1), on the Nishimori line.

rho = 0.5;
alphaMin = 0.4;
alphaMax = 0.6;
dAlpha = 0.02;
accuracyE = 1e-15;
dump = 0.5;
varNoise = 0.01;

fa = @(S2, R, rho) tanh(R ./ S2);
% fa = @(S2, R, rho_) (rho_ .* exp(-(1 - R).^2 ./ (2 .* S2) ) - (1 - rho_) .* exp(-(1 + R).^2 ./ (2 * S2) ) ) ./ (rho_ .* exp(-(1 - R).^2 ./ (2 .* S2) ) + (1 - rho_) .* exp(-(1 + R).^2 ./ (2 .* S2) ) );

% si p(x) = rho * delta(x - 1) + (1 - rho) * delta(x)
% fa = @(S2, R, rho) rho ./ (rho + (1 - rho) .* exp((-R.^2 + (1 - R).^2) ./ (2 * S2) ) );

G = @(z) exp(-z.^2 / 2) / sqrt(2 * pi);

alpha = alphaMin; i = 1;
while (alpha <= alphaMax)
    
    E = 1; dE = accuracyE + 1;
    while (dE >= accuracyE)
        
        int1 = integral(@(z) G(z) .* (fa((varNoise + E) ./ alpha, 1 + z .* sqrt((E + varNoise) ./ alpha), rho) - 1).^2, -10, 10);
        int2 = integral(@(z) G(z) .* (fa((varNoise + E) ./ alpha, -1 + z .* sqrt((E + varNoise) ./ alpha), rho) + 1).^2, -10, 10);
        Enew = rho * int1 + (1 - rho) * int2;
        
        % si p(x) = rho * delta(x - 1) + (1 - rho) * delta(x)
        %                 int1 = integral(@(z) G(z) .* fa((varNoise + E) ./ alpha, z .* sqrt((E + varNoise) ./ alpha), rho).^2, -10, 10);
        %                 int2 = integral(@(z) G(z) .* (fa((varNoise + E) ./ alpha, 1 + z .* sqrt((E + varNoise) ./ alpha), rho) - 1).^2, -10, 10);
        %                 Enew = rho * int2 + (1 - rho) * int1;
        
        dE = abs(Enew - E);
        E = dump * E + (1 - dump) * Enew;
    end
    MseTab(i) = E;
    
%     E = 1; dE = accuracyE + 1;
%     while (dE >= accuracyE)
%         
%         Enew = 2 .* (1 - erf(sqrt(alpha ./ (2 .* (varNoise + E) ) ) ) );
%         
%         dE = abs(Enew - E);
%         E = dump * E + (1 - dump) * Enew;
%     end
    MseSignTab(i) = 2 .* (1 - erf(sqrt(alpha ./ (2 .* (varNoise + E) ) ) ) );
    
    alphaTab(i) = alpha;
    alpha = alpha + dAlpha; i = i + 1;
end

plot(alphaTab, MseTab,'b*',alphaTab, MseSignTab,'r*')