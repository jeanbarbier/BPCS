function [F_E_result] = F_E(alpha,rho,var_noise,var_psi,E_min,E_max,delta_E)

E=E_min;
i=1;

while (E <= E_max)    
    f=( var_psi.*( var_psi - 1 + ( var_noise + E )./alpha ) - ( var_noise + E)./alpha )./( ( var_psi + ( var_noise + E )./alpha ).*( 1 + ( var_noise + E )./alpha ) );    
    g=( var_psi.*( 1 + ( var_noise + E )./alpha ) - ( var_noise + E )./alpha - 1 )./( ( var_psi + ( var_noise + E )./alpha ).*( 1 + ( var_noise + E )./alpha ) );
        
    integrand_1=@(u) exp(-(u.^2)./2).*((2.*pi).^(-1./2)).*log( (1-rho).*( (var_noise + E)./(var_noise + E + var_psi.*alpha) ).^(1./2).*exp(0.5.*(u.^2).*f) + rho.*( (var_noise + E)./(var_noise + E + alpha) ).^(1./2) );    
    integrand_2=@(u) exp(-(u.^2)./2).*((2.*pi).^(-1./2)).*log( (1-rho).*( (var_noise + E)./(var_noise + E + var_psi.*alpha) ).^(1./2).*exp(0.5.*(u.^2).*g) + rho.*( (var_noise + E)./(var_noise + E + alpha) ).^(1./2) );    
    
    [int_1,err_1]=quadgk(integrand_1,-10,10);   
    [int_2,err_2]=quadgk(integrand_2,-10,10);   
    
    part_0=-alpha.*0.5.*( 1 + log(var_noise + E) + (rho + (1-rho).*var_psi - E)./(var_noise + E) );    
    part_1=(1-rho).*( 0.5.*(1 + (var_psi.*alpha)./(var_noise + E))./(1 + (var_noise + E)./alpha) + int_1 );    
    part_2=rho.*( 0.5.*(1 + alpha./(var_noise + E))./(1 + (var_noise + E)./alpha) + int_2 );   
    
    F_E_result(1,i)=part_0 + part_1 + part_2;
    E_plot(1,i)=E;    
    i=i+1;
    E=E.*delta_E;   
end

semilogx(E_plot,F_E_result);
    
end

