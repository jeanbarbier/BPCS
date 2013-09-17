classdef Dump
    
    properties
        f; val; func;
    end
    
    methods
        
        function dump = Dump(dump_init,choice_prior)
            dump.f = []; dump.val = dump_init; 
            switch choice_prior
                case 'SparseGauss'
                    dump.func = 'dumpSG';
                case 'SparseGaussCut'
                    dump.func = 'dumpSG';
                case '2Gauss'
                    dump.func = 'dump2G';
                case 'SparseExponential'
                    dump.func = 'dumpSE';
                case 'SparseConstant'
                    dump.func = 'dumpSC';
                case 'SparseBinary'
                    dump.func = 'dumpSB';
                case 'L1'
                    dump.func = 'dumpL1';
                case 'Laplace'
                    dump.func = 'dumpLap';
                case 'SGMF'
                    dump.func = 'dumpSGMF';
            end
        end
        
        function dump = dumpSG(dump,rho,M,N,R,S2,G,G2,var_noise,V,Y,W,f_a,f_c,varargin)
            m = varargin{1}; var_gauss = varargin{2};           
            log_Z_mu = -0.5 .* ( (Y - W).^2 ./ (var_noise + V) + log(2 .* pi .* (var_noise + V) ) );
            log_Z_i1 = R.^2 ./ (2 .* S2) + log( (1 - rho) .* exp(-0.5 .* R.^2 ./ S2) + rho .* sqrt(S2 ./ (S2 + var_gauss) ) .* exp(-0.5 .* (R - m).^2 ./ (S2 + var_gauss) ) );
            log_Z_i2 = -f_a .* R ./ S2 + 0.5 ./ S2 .* (f_c + f_a.^2) + 0.5 .* f_c .* (((Y - W).^2 .* (var_noise + V).^(-2) ) * G2);    
            dump.f = (sum(log_Z_mu) + sum(log_Z_i1 + log_Z_i2) ) ./ N;
        end
        
        function dump = dumpSGMF(dump,rho,M,N,R,S2,G,G2,var_noise,V,Y,W,f_a,f_c,varargin)
            m = varargin{1}; var_gauss = varargin{2}; 
            Gfa = G * f_a';
            part1 = -0.5 .* sum(Y.^2 - 2 .* Y .* Gfa' + (Gfa.^2)' + (G2 * f_c')' );
            part2 = var_noise .* ((f_c + f_a.^2) * ((2 .* S2).^(-1))' - f_a * (R ./ S2)');
            part3 = var_noise .* sum(log((1 - rho) .* exp(-0.5 ./ (S2 .* var_gauss) .* (m .* S2 + R .* var_gauss).^2 ./ (S2 + var_gauss) ) + rho .* exp(-0.5 .* m.^2 ./ var_gauss) .* sqrt(S2 ./ (S2 + var_gauss) ) ) );
            part4 = var_noise .* 0.5 ./ (var_gauss .* S2) * ((m .* S2 + R .* var_gauss).^2 ./ (S2 + var_gauss) )';
            dump.f = (part1 + part2 + part3 + part4) ./ N;
        end
        
        function dump = dumpSB(dump,rho,M,N,R,S2,G,G2m,var_noise,V,Y,W,f_a,f_c,varargin)
            a = exp((-0.5 + R) ./ S2);
            b = a .* (R ./ S2 - 0.5 .* (1 ./ S2 + G2m .* sum((Y - W).^2 .* (var_noise + V).^(-2)) ) );
            log_Z_mu = -0.5 .* ( (Y - W).^2 ./ (var_noise + V) + log(2 .* pi .* (var_noise + V) ) );
            log_Z_i1 = log( (1 - rho) + rho .* a);
            log_Z_i2 = -rho  .*  b ./ exp(log_Z_i1);
            dump.f = sum(log_Z_mu) + sum(log_Z_i1 + log_Z_i2) ./ N;
        end
        
    end
    
end