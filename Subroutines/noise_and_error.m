classdef noise_and_error
    % This class contains the noise learning and error estimations
    
    properties
        true_error; convergence; var_noise_old; var_noise; conv_crit; dump_learn;
    end
    
    methods
        
        function obj = noise_and_error(conv_crit,var_noise_init,dump_learn)
            % Constructor function
            obj.conv_crit = conv_crit; obj.var_noise_old = var_noise_init; obj.var_noise = var_noise_init; obj.dump_learn = dump_learn;
        end
        
        function obj = learn_noise(obj,Y,W,V)
            obj.var_noise = obj.dump_learn .* obj.var_noise_old * (1 - obj.dump_learn) .* ((Y - W).^2 * (1 + V ./ obj.var_noise_old).^(-2)' ) ./ sum((1 + V ./ obj.var_noise_old).^(-1) );
            if (obj.var_noise < 1e-100); obj.var_noise = 1e-100; end;
            obj.var_noise_old = obj.var_noise;
        end
        
        function obj = compute_true_MSE(obj,signal,X)
            obj.true_error = sum((signal - X).^2) ./ max(size(X) );
        end
        
        function obj = compute_convergence(obj,X_old,X)
            obj.convergence = sum(vabs(X_old - X) ) ./ max(size(X) );
        end
        
    end
    
end