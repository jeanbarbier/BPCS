classdef Prior
    % This class contains all the prior-dependent functions including learnings
    
    properties
        av_mess; av_mess_old; var_mess; var_mess_old; R; S2; rho; learn; N; alpha; func; dump_learn; t; method;
        param_1; param_2; param_3; param_4;
        %        param_1 = m_gauss; param_2 = var_gauss; % Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / (sqrt(2 * pi) * dev) * exp(-(x - m)^2 / (2 * dev^2) )
        %        param_1 = m_gauss; param_2 = var_gauss; param_3 = G; p(x) ~ [(1 - rho) * delta(x) + rho / sqrt(2 * pi * var) * exp(-(x - m)^2 / (2 * var) )] * I(|x| < G)
        %        param_1 = m_1; param_2 = m_2; param_3 = var_1; param_4 = var_2; % Mixture of two gaussians : p(x) ~ (1 - rho) * / (sqrt(2 * pi * var_1) ) * exp(-(x - m_1)^2 / (2 * var_1) ) + rho / (sqrt(2 * pi * var_2) ) * exp(-(x - m_2)^2 / (2 * var_2) )
        %        param_1 = expo; % Exponential sparse prior : p(x) ~ (1 - rho) * delta(x) + rho * I(x > 0) * exp(-expo * x), expo > 0
        %        param_1 = c_down; param_2 = c_up; % Unity inside a finite interval sparse prior : p(x) ~ (1 - rho) * delta(x) + rho * I(c_down < x < c_up)
        %        param_1 = beta % Laplace prior : p(x) ~ 2 / beta * exp{-beta * |x|}
        %        param_1 = x_mean; param_2 = var_noise; % Complex prior
    end
    
    methods
        
        function prior = Prior(rho_init, N, alpha, learn, choice_prior, dump_learn, R_init, S2_init, av_mess_init, var_mess_init, method, varargin)
            % Constructor function
            prior.R = R_init; prior.S2 = S2_init; prior.rho = rho_init; prior.learn = learn; prior.N = N; prior.alpha = alpha; prior.av_mess = av_mess_init; prior.av_mess_old = av_mess_init; prior.var_mess = var_mess_init; prior.var_mess_old = var_mess_init; prior.dump_learn = dump_learn; prior.method = method;
            switch choice_prior
                case 'SparseGauss'
                    prior.param_1 = varargin{1};
                    prior.param_2 = varargin{2};
                    prior.func = 'PriorSG';
                    disp('SparseGauss')
                case 'SparseGaussCut'
                    prior.param_1 = varargin{1};
                    prior.param_2 = varargin{2};
                    prior.param_3 = varargin{3};
                    prior.func = 'PriorSGC';
                    disp('SparseGaussCut')
                case 'SparseGaussPositive'
                    prior.param_1 = varargin{1};
                    prior.param_2 = varargin{2};
                    prior.param_3 = varargin{3};
                    prior.func = 'PriorSGP';
                    disp('SparseGaussPositive')
                case '2Gauss'
                    prior.param_1 = varargin{1};
                    prior.param_2 = varargin{2};
                    prior.param_3 = varargin{3};
                    prior.param_4 = varargin{4};
                    prior.func = 'Prior2G';
                    disp('2Gauss')
                case 'SparseExponential'
                    prior.param_1 = varargin{1};
                    prior.func = 'PriorSE';
                    disp('SparseExponential')
                case 'SparseConstant'
                    prior.param_1 = min(varargin{1},varargin{2});
                    prior.param_2 = max(varargin{1},varargin{2});
                    prior.func = 'PriorSC';
                    disp('SparseConstant')
                case 'SparseBinary'
                    prior.func = 'PriorSB';
                    disp('SparseBinary')
                case 'L1'
                    prior.param_1 = min(varargin{1},varargin{2});
                    prior.param_2 = max(varargin{1},varargin{2});
                    prior.func = 'PriorL1';
                    disp('L1')
                case 'Laplace'
                    prior.param_1 = varargin{1};
                    prior.func = 'PriorLap';
                    disp('Laplace')
                case 'Binary1'
                    prior.func = 'PriorPm1';
                    disp('PriorPm1')
                case 'Complex'
                    prior.param_1 = varargin{1};
                    prior.param_2 = varargin{2};
                    prior.func = 'PriorComplex';
                    disp('Complex')
                otherwise
                    disp('unknown prior')
            end
        end
        
        function prior = PriorSG(prior)
            % Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var) * exp(-(x - m)^2 / (2 * var) )
            R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho; m_ = prior.param_1; var_ = prior.param_2; N_ = prior.N;
            prior.av_mess_old = prior.av_mess;
            a = exp(-R_.^2 ./ (2 .* S2_) + 0.5 .* (m_ - R_).^2 ./ (var_ + S2_) );
            c = 1 ./ sqrt(2 .* pi .* (S2_ + var_) );
            Z = (1 - rho_) ./ sqrt(2 .* pi .* S2_) .* a + rho_ .* c;
            prior.av_mess = (rho_ .* c .* (m_ .* S2_ + R_ .* var_) ./ (S2_ + var_) ) ./ Z;
            f_b = (rho_ ./ sqrt(2 .* pi .* S2_ .* var_) .* (S2_.^(-1) + var_.^(-1) ).^(-3 ./ 2) .* (1 + (m_ ./ var_ + R_ ./ S2_).^2 ./ (1 ./ S2_ + 1 ./ var_) ) ) ./ Z;
            prior.var_mess = max(1e-18,f_b - prior.av_mess.^2);
            if (prior.learn == 1)
                a = exp(-R_.^2 ./ (2 .* S2_) );
                Z_rho = (1 - rho_) .* a ./ sqrt(2 .* pi .* S2_) + prior.av_mess .* (var_ + S2_) ./ (m_ .* S2_ + R_ .* var_);
                prior.rho = prior.dump_learn .* prior.rho + (1 - prior.dump_learn) .* abs(((prior.av_mess .* (var_ + S2_) ./ (m_ .* S2_ + R_ .* var_) ) * Z_rho.^(-1)') ./ (a ./ sqrt(2 .* pi .* S2_) * Z_rho.^(-1)') );
                prior.param_1 = prior.dump_learn .* prior.param_1 + (1 - prior.dump_learn) .* sum(prior.av_mess) ./ (rho_ .* N_);
                prior.param_2 = prior.dump_learn .* prior.param_2 + (1 - prior.dump_learn) .* max(1e-17,sum(f_b) ./ (rho_ .* N_) - m_.^2);
                if (prior.rho > prior.alpha); prior.rho = prior.alpha; end;
            end
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = PriorSGC(prior)
            % SparseGauss prior with a cut on the values of x : p(x) ~ [(1 - rho) * delta(x) + rho / sqrt(2 * pi * var) * exp(-(x - m)^2 / (2 * var) )] * I(|x| < G)
            R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho; m_ = prior.param_1; var_ = prior.param_2; G_ = prior.param_3;
            prior.av_mess_old = prior.av_mess;
            Vp = erf(((G_ + R_) .* var_ + (G_ + m_) .* S2_) ./ sqrt(2 .* S2_ .* var_ .* (S2_ + var_) ) );
            Vm = erf(((-G_ + R_) .* var_ + (-G_ + m_) .* S2_) ./ sqrt(2 .* S2_ .* var_ .* (S2_ + var_) ) );
            Kp = sqrt(var_ .* S2_) ./ (var_ + S2_) .* exp(-0.5 .* ((R_ + G_).^2 ./ S2_ + (m_ + G_).^2 ./ var_) );
            Km = sqrt(var_ .* S2_) ./ (var_ + S2_) .* exp(-0.5 .* ((R_ - G_).^2 ./ S2_ + (m_ - G_).^2 ./ var_) );
            Fp = sqrt(pi ./ 2) .* (R_ .* var_ + m_ .* S2_) ./ (S2_ + var_).^(3 ./ 2) .* erf((R_ .* var_ + m_ .* S2_ + G_ .* (var_ + S2_) ) ./ sqrt(2 .* var_ .* S2_ .* (S2_ + var_) ) ) .* exp(-0.5 .* (m_ - R_).^2 ./ (S2_ + var_) );
            Fm = sqrt(pi ./ 2) .* (R_ .* var_ + m_ .* S2_) ./ (S2_ + var_).^(3 ./ 2) .* erf((R_ .* var_ + m_ .* S2_ - G_ .* (var_ + S2_) ) ./ sqrt(2 .* var_ .* S2_ .* (S2_ + var_) ) ) .* exp(-0.5 .* (m_ - R_).^2 ./ (S2_ + var_) );
            Ep = erf((-R_ .* var_ - m_ .* S2_ + G_ .* (var_ + S2_) ) ./ sqrt(2 .* var_ .* S2_ .* (S2_ + var_) ) );
            Em = erf((-R_ .* var_ - m_ .* S2_ - G_ .* (var_ + S2_) ) ./ sqrt(2 .* var_ .* S2_ .* (S2_ + var_) ) );
            Gp = -sqrt(S2_ .* var_) ./ (S2_ + var_).^2 .* (R_ .* var_ + m_ .* S2_ + G_ .* (S2_ + var_) ) .* exp(-0.5 .* (R_.^2 ./ S2_ + m_.^2 ./ var_ + G_.^2 .* (S2_ + var_) ./ (S2_ .* var_) - 2 .* G_ .* m_ ./ var_ - 2 .* G_ .* R_ ./ S2_) );
            Gm = -sqrt(S2_ .* var_) ./ (S2_ + var_).^2 .* (R_ .* var_ + m_ .* S2_ - G_ .* (S2_ + var_) ) .* exp(-0.5 .* (R_.^2 ./ S2_ + m_.^2 ./ var_ + G_.^2 .* (S2_ + var_) ./ (S2_ .* var_) + 2 .* G_ .* m_ ./ var_ + 2 .* G_ .* R_ ./ S2_) );
            U = sqrt(pi ./ 2) .* ((R_ .* var_ + m_ .* S2_).^2 + S2_.^2 .* var_ + var_.^2 .* S2_) ./ (S2_ + var_).^(5 ./ 2) .* exp(-0.5 .* (m_ - R_).^2 ./ (S2_ + var_) );
            Z = (1 - rho_) .* exp(-0.5 .* R_.^2 ./ S2_) ./ sqrt(2 .* pi .* S2_) + rho_ .* exp(-0.5 .* (m_ - R_).^2 ./ (S2_ + var_) ) ./ sqrt(8 .* pi .* (S2_ + var_) ) .* (Vp - Vm);
            prior.av_mess = rho_ ./ (Z .* 2 .* pi) .* (Kp + Fp - Km - Fm);
            f_b = rho_ ./ (Z .* 2 .* pi) .* (Gp - Gm + U .* (Ep - Em) );
            prior.var_mess = max(1e-18,f_b - prior.av_mess.^2);
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = PriorSGP(prior)
            % Positive Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var) * exp(-(x - m)^2 / (2 * var) ) * I(x > 0)
            R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho; m_ = prior.param_1; var_ = prior.param_2; N_ = prior.N;
            prior.av_mess_old = prior.av_mess;
            aa = -(m_ .* S2_ + R_ .* var_) ./ sqrt(2 .* S2_ .* var_ .* (S2_ + var_) );
            A = erfc(aa);
            B = exp(-1 ./ (2 .* var_ .* S2_) .* (m_.^2 .* S2_ + R_.^2 .* var_ - (m_ .* S2_ + R_ .* var_).^2 ./ (S2_ + var_) ) );
            C = exp(-(m_.^2 .* S2_ + R_.^2 .* var_) ./ (2 .* var_ .* S2_) );
            Z = (1 - rho_) .* exp(-R_.^2 ./ (2 .* S2_) ) + rho_ .* sqrt(S2_ .* var_ .* pi .* .5 ./ (S2_ + var_) ) .* A .* B;
            first = rho_ .* S2_ .* var_ ./ (S2_ + var_) .* C + rho_ .* sqrt(S2_ .* var_ .* pi .* .5 ./ (S2_ + var_).^3) .* B .* A .* (m_ .* S2_ + R_ .* var_);
            second = rho_ .* B .* sqrt(S2_ .* var_ .* pi .* .5 ./ (S2_ + var_).^3) .* A .* (S2_ .* var_ + (m_ .* S2_ + R_ .* var_).^2 ./ (S2_ + var_) ) + rho_ .* C .* S2_ .* var_ .* (m_ .* S2_ + R_ .* var_) ./ (S2_ + var_).^2;
            prior.av_mess = first ./ Z;
            prior.var_mess = max(1e-18,second ./ Z - prior.av_mess.^2);
            if (prior.learn == 1)
                m_up = (sqrt(S2_ ./ (S2_ + var_) ) .* B .* (A .* (m_ .* S2_ + R_ .* var_) ./ (var_ .* (S2_ +  var_) ) + sqrt(2 .* S2_ ./ (pi .* var_ .* (S2_ + var_) ) ) .* exp(-aa.^2) ) ) * Z.^(-1)';
                m_down = (sqrt(S2_ ./ (S2_ + var_) ) .* B .* A ./ var_) * Z.^(-1)';
                %                 rho_down = (sqrt(S2_ .* var_ .* pi .* .5 ./ (S2_ + var_) ) .* A .* B) * Z.^(-1)';
                %                 rho_up = exp(-.5 .* R_.^2 ./ S2_) * Z.^(-1)';
                %                 prior.rho = prior.dump_learn .* prior.rho + (1 - prior.dump_learn) .* rho_ .* rho_up ./ rho_down;
                %                 prior.param_1 = prior.dump_learn .* prior.param_1 + (1 - prior.dump_learn) .* sum(prior.av_mess) ./ (rho_ .* N_);
                prior.param_1 = prior.dump_learn .* prior.param_1 + (1 - prior.dump_learn) .* m_up ./ m_down;
                prior.param_2 = prior.dump_learn .* prior.param_2 + (1 - prior.dump_learn) .* max(1e-18,sum(second ./ Z) ./ (rho_ .* N_) - m_.^2);
                if (prior.rho > prior.alpha); prior.rho = prior.alpha; end;
            end
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = Prior2G(prior)
            % Mixture of two gaussians : p(x) ~ (1 - rho) * / (sqrt(2 * pi * var_1) ) * exp(-(x - m_1)^2 / (2 * var_1) ) + rho / (sqrt(2 * pi * var_2) ) * exp(-(x - m_2)^2 / (2 * var_2) )
            R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho; m_1_ = prior.param_1; m_2_ = prior.param_2; var_1_ = prior.param_3; var_2_ = prior.param_4; N_ = prior.N;
            prior.av_mess_old = prior.av_mess;
            %             a = exp(-0.5 .* (m_2_ - R_).^2 ./ (var_2_ + S2_) );
            %             b = exp(-0.5 .* (R_ - m_1_).^2 ./ (var_1_ + S2_) );
            %             c = 1 ./ sqrt(2 .* pi .* (S2_ + var_1_) );
            %             d = 1 ./ sqrt(2 .* pi .* (S2_ + var_2_) );
            %             e = (m_2_ .* S2_ + R_ .* var_2_) ./ (S2_ + var_2_);
            %             f = (R_ .* var_1_ + m_1_ .* S2_) ./ (S2_ + var_1_);
            %             Z = rho_ .* a .* d + (1 - rho_) .* b .* c;
            %             f_b = ((1 - rho_) ./ sqrt(2 .* pi .* S2_ .* var_1_) .* (1 ./ var_1_ + 1 ./ S2_).^(-3 ./ 2) .* b .* (1 + (R_ ./ S2_).^2 ./ (1 ./ var_1_ + 1 ./ S2_) ) + rho_ ./ sqrt(2 .* pi .* S2_ .* var_2_) .* (1 ./ var_2_ + 1 ./ S2_).^(-3 ./ 2) .* a .* (1 + (m_2_ ./ var_2_ + R_ ./ S2_).^2 ./ (1 ./ var_2_ + 1 ./ S2_) ) ) ./ Z;
            %             prior.av_mess = (rho_ .* a .* d .* e + (1 - rho_) .* b .* c .* f) ./ Z;
            %             prior.var_mess = max(1e-18,f_b - prior.av_mess.^2);
            
            % Another equivalent form
            expo1 = exp(-.5 .* (m_1_ - R_).^2 ./ (var_1_ + S2_) );
            expo2 = exp(-.5 .* (m_2_ - R_).^2 ./ (var_2_ + S2_) );
            f1 = 1 ./ sqrt(2 .* pi .* (var_1_ + S2_) ) .* expo1;
            f2 = 1 ./ sqrt(2 .* pi .* (var_2_ + S2_) ) .* expo2;
            g1 = (m_1_ .* S2_ + R_ .* var_1_) ./ (var_1_ + S2_) .* f1;
            g2 = (m_2_ .* S2_ + R_ .* var_2_) ./ (var_2_ + S2_) .* f2;
            h1 = 1 ./ sqrt(2 .* pi .* var_1_ .* S2_) .* (1 ./ var_1_ + 1 ./ S2_).^(-3 ./ 2) .* (1 + (m_1_ ./ var_1_ + R_ ./ S2_).^2 ./ (1 ./ var_1_ + 1 ./ S2_) ) .* expo1;
            h2 = 1 ./ sqrt(2 .* pi .* var_2_ .* S2_) .* (1 ./ var_2_ + 1 ./ S2_).^(-3 ./ 2) .* (1 + (m_2_ ./ var_2_ + R_ ./ S2_).^2 ./ (1 ./ var_2_ + 1 ./ S2_) ) .* expo2;
            Z = (1 - rho_) .* f1 + rho_ .* f2;
            f_b = ((1 - rho_) .* h1 + rho_ .* h2) ./ Z;
            prior.av_mess = ((1 - rho_) .* g1 + rho_ .* g2) ./ Z;
            prior.var_mess = max(1e-18, f_b - prior.av_mess.^2);
            if (prior.learn == 1)
                %                 U = 1 ./ S2_;
                %                 V = R_ ./ S2_;
                %                 tho = (1 - rho_) .* (var_1_ .* U + 1).^(-1 ./ 2) .* exp(V.^2 ./ (2 .* (U + 1 ./ var_1_) ) - ((V + m_2_ ./ var_2_).^2) ./ (2 .* (U + 1 ./ var_2_) ) ) + rho_ .* (var_2_ .* U + 1).^(-1 ./ 2) .* exp(-m_2_.^2 ./ (2 .* var_2_) );
                %                 up = (rho_ .* exp(-m_2_.^2 ./ (2 .* var_2_) ) ./ sqrt(var_2_ .* U + 1) ) * tho.^(-1)';
                %                 down = (exp(V.^2 ./ (2 .* (U + 1 ./ var_1_) ) - ((V + m_2_ ./ var_2_).^2) ./ (2 .* (U + 1 ./ var_2_) ) ) ./ sqrt(var_1_ .* U + 1) ) * tho.^(-1)';
                %                 prior.rho = prior.dump_learn .* prior.rho + (1 - prior.dump_learn) .* up ./ down;
                %                 if (prior.rho > prior.alpha); prior.rho = prior.alpha; end;
                %                 up_var_2 = (U + 1 ./ var_2_).^(-1 ./ 2) * ((m_2_.^2 - 2 .* m_2_ .* ((V + m_2_ ./ var_2_) ./ (U + 1 ./ var_2_) ) + (((V + m_2_ ./ var_2_) ./ (U + 1 ./ var_2_) ).^2) + (U + 1 ./ var_2_).^(-1) ) .* tho.^(-1) )';
                %                 down_var_2 = (U + 1 ./ var_2_).^(-1 ./ 2) * tho.^(-1)';
                %                 prior.param_3 = prior.dump_learn .* prior.param_3 + (1 - prior.dump_learn) .* up_var_2 ./ down_var_2;
                %                 up_var_1 = ((U + 1 ./ var_1_).^(-1 ./ 2) .* exp(V.^2 ./ (2 .* (U + 1 ./ var_1_) ) - ((V + m_2_ ./ var_2_).^2) ./ (2 .* (U + 1 ./ var_2_) ) ) .* (((V ./ (U + 1 ./ var_1_) ).^2) + (U + 1 ./ var_1_).^(-1) ) ) * tho.^(-1)';
                %                 down_var_1 = ((U + 1 ./ var_1_).^(-1 ./ 2) .* exp(V.^2 ./ (2 .* (U + 1 ./ var_1_) ) - ((V + m_2_ ./ var_2_).^2) ./ (2 .* (U + 1 ./ var_2_) ) ) ) * (tho.^(-1))';
                %                 prior.param_2 = prior.dump_learn .* prior.param_2 + (1 - prior.dump_learn) .* up_var_1 ./ down_var_1;
                %                 prior.param_1 = prior.dump_learn .* prior.param_1 + (1 - prior.dump_learn) .* 1 ./ (N_ .* rho_) .* sum(prior.av_mess);
                %             end
                
                % Another equivalent form
                prior.rho = prior.dump_learn .* rho_ + (1 - prior.dump_learn) .* rho_ .* (f2 * Z.^(-1)') ./ (f1 * Z.^(-1)');
                if (prior.rho > prior.alpha); prior.rho = prior.alpha; end;
                %                 prior.param_1 = prior.dump_learn .* m_1_ + (1 - prior.dump_learn) .* (mean(prior.av_mess) - rho_ .* m_2_) ./ (1 - rho_);
                %                 prior.param_2 = prior.dump_learn .* m_2_ + (1 - prior.dump_learn) .* (mean(prior.av_mess) - (1 - rho_) .* m_1_) ./ rho_;
                %                 prior.param_1 = prior.dump_learn .* m_1_ + (1 - prior.dump_learn) .* ((f1 .* R_ ./ (var_1_ + S2_) ) * Z.^(-1)') ./ ((f1 ./ (var_1_ + S2_) ) * Z.^(-1)');
                prior.param_2 = prior.dump_learn .* m_2_ + (1 - prior.dump_learn) .* ((f2 .* R_ ./ (var_2_ + S2_) ) * Z.^(-1)') ./ ((f2 ./ (var_2_ + S2_) ) * Z.^(-1)');
                prior.param_3 = prior.dump_learn .* var_1_ + (1 - prior.dump_learn) .* var_1_ .* ((expo1 .* (m_1_ - R_).^2 ./ (var_1_ + S2_).^(5 ./ 2) ) * Z.^(-1)') .* (expo1 ./ (var_1_ + S2_).^(3 ./ 2) * Z.^(-1)').^(-1);
                prior.param_4 = prior.dump_learn .* var_2_ + (1 - prior.dump_learn) .* var_2_ .* ((expo2 .* (m_2_ - R_).^2 ./ (var_2_ + S2_).^(5 ./ 2)) * Z.^(-1)') .* (expo2 ./ (var_2_ + S2_).^(3 ./ 2) * Z.^(-1)').^(-1);
            end
            %             if (strcmp(prior.method,'AMPhadamard') & (prior.t < 100) )
            %                 prior.av_mess(1) = prior.av_mess_old(1);
            %                 prior.var_mess(1) = prior.var_mess_old(1);
            %             end
        end
        
        function prior = PriorSB(prior)
            % Binary prior
            R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho; N_ = prior.N;
            prior.av_mess_old = prior.av_mess;
            prior.av_mess = rho_ ./ (rho_ + (1 - rho_) .* exp(  (1 - 2 .* R_) ./ (2 .* S2_) ) );
            prior.var_mess = prior.av_mess .* (1 - prior.av_mess);
            if (prior.learn == 1)
                prior.rho = prior.dump_learn .* prior.rho + (1 - prior.dump_learn) .* sum(prior.av_mess) ./ N_;
                if (prior.rho > prior.alpha); prior.rho = prior.alpha; end;
            end
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = PriorSE(prior)
            % Exponential sparse prior : p(x) ~ (1 - rho) * delta(x) + rho * I(x > 0) * exp(-expo * x), expo > 0
            R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho; expo_ = prior.param_1; N_ = prior.N;
            prior.av_mess_old = prior.av_mess;
            a = exp(-R_.^2 ./ (2 .* S2_) );
            b = exp(-expo_ .* R_ + expo_.^2 .* S2_ ./ 2);
            c = erfc((expo_ .* sqrt(S2_) - R_ ./ sqrt(S2_) ) ./ sqrt(2) );
            Z = (1 - rho_) ./ sqrt(2 .* pi .* S2_) .* a + rho_ .* b ./ 2 .* c;
            prior.av_mess = rho_ .*  (sqrt(S2_ ./ (2 .* pi) ) .* a + (R_ - expo_ .* S2_) ./ 2 .* b .* c) ./ Z;
            f_b = rho_ .* (S2_ ./ sqrt(2 .* pi) .* (-expo_ .* sqrt(S2_) + R_ ./ sqrt(S2_) ) .* a + b ./ 2 .* c .* (S2_ + (R_ - expo_ .* S2_).^2 ) ) ./ Z;
            prior.var_mess = max(1e-18,f_b - prior.av_mess.^2);
            if (prior.learn == 1)
                Z_rho = (1 - rho_) .* a ./ sqrt(2 .* pi .* S2_) + (prior.av_mess - rho_ .* sqrt(S2_ ./ (2 .* pi) .* a) ) ./ (R_ - expo_ .* S2_);
                prior.rho = prior.dump_learn .* prior.rho + (1 - prior.dump_learn) .* abs(((prior.av_mess - rho_ .* sqrt(S2_ ./ (2 .* pi) .* a) ) ./ (R_ - expo_ .* S2_) * Z_rho.^(-1)') ./ (a ./ sqrt(2 .* pi .* S2_) * Z_rho.^(-1)') );
                prior.param_1 = prior.dump_learn .* prior.param_1 + (1 - prior.dump_learn) .* sqrt(rho_ .* N_ ./ sum(prior.av_mess) );
                if (prior.rho > prior.alpha); prior.rho = prior.alpha; end;
            end
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = PriorSC(prior)
            % Unity inside a finite interval sparse prior : p(x) ~ (1 - rho) * delta(x) + rho * I(c_down < x < c_up)
            R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho; c_down_ = prior.param_1; c_up_ = prior.param_2;
            prior.av_mess_old = prior.av_mess;
            a = exp(-R_.^2 ./ (2 .* S2_) );
            b = erfc((R_ - c_down_) ./ sqrt(S2_ .* 2) );
            c = erfc((R_ - c_up_) ./ sqrt(S2_ .* 2) );
            d = exp(-0.5 .* (c_down_ - R_).^2 ./ S2_ );
            e = exp(-0.5 .* (c_up_ - R_).^2 ./ S2_);
            f = 1 ./ R_ .* (prior.av_mess - rho_ .* S2_ ./ sqrt(2 .* pi) .* (d - e) );
            Z = (1 - rho_) ./ sqrt(2 .* pi .* S2_) .* a + rho_ .* 0.5 .* (c - b);
            prior.av_mess = rho_ .* (sqrt(S2_ ./ (2 .* pi) ) .* (d - e) + R_ ./ 2 .* (c - b) ) ./ Z;
            f_b = rho_ .* (sqrt(S2_ ./ pi) .* (d .* (1 ./ sqrt(2) .* (c_down_ - R_) + sqrt(2) .* R_) - e .* (1 ./ sqrt(2) .* (c_up_ - R_) + sqrt(2) .* R_) ) + 0.5 .* (R_.^2 + S2_) .* (c - b) ) ./ Z;
            prior.var_mess = max(1e-18,f_b - prior.av_mess.^2);
            if (prior.learn == 1)
                Z_rho = (1 - rho_) .* a ./ sqrt(2 .* pi .* S2_) + f;
                prior.rho = prior.dump_learn .* prior.rho + (1 - prior.dump_learn) .* abs((f * Z_rho.^(-1)') ./ (a ./ sqrt(2 .* pi .* S2_) * Z_rho.^(-1)') );
                if (prior.rho > prior.alpha); prior.rho = prior.alpha; end;
            end
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = PriorL1(prior)
            % L1 optimization : p(x) ~ lim_{beta -> infinity} exp{-beta * |x|}
            prior.av_mess_old = prior.av_mess; min_ = prior.param_1; max_ = prior.param_2;
            R_ = prior.R; S2_ = prior.S2;
            prior.av_mess = min(max_,(R_ > 0) .* (R_ - S2_) .* (R_ > S2_) ) + max(min_,(R_ < 0) .* (R_ + S2_) .* (-R_ > S2_) );
            prior.var_mess = S2_ .* (abs(R_) > S2_);
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = PriorLap(prior)
            prior.av_mess_old = prior.av_mess; beta_ = prior.param_1;
            R_ = prior.R; S2_ = prior.S2;
            erfc_p = erfc((R_ + beta_ .* S2_) ./ sqrt(2 .* S2_) );
            erfc_m = erfc((-R_ + beta_ .* S2_) ./ sqrt(2 .* S2_) );
            z = erfc_p + erfc_m .* exp(-2 .* beta_ .* R_);
            f_b_part1 = -4 .* beta_ .* S2_.^(3 ./ 2) ./ (sqrt(2 .* pi) .* (exp((beta_ .* S2_ + R_).^2 ./ (2 .* S2_) ) .* erfc_p + exp((beta_ .* S2_ - R_).^2 ./ (2 .* S2_) ) .* erfc_m) );
            f_b_part2 = (((R_ + beta_ .* S2_).^2 + S2_) .* erfc_p + ((R_ - beta_ .* S2_).^2 + S2_) .* erfc_m .* exp(-2 .* beta_ .* R_) ) ./ z;
            prior.av_mess = ((R_ + beta_ .* S2_) .* erfc_p + (R_ - beta_ .* S2_) .* erfc_m .* exp(-2 .* beta_ .* R_) ) ./ z;
            prior.var_mess = max(1e-18,f_b_part1 + f_b_part2 - prior.av_mess.^2);
            if (strcmp(prior.method,'AMPhadamard') )
                prior.av_mess(1) = prior.av_mess_old(1);
                prior.var_mess(1) = prior.var_mess_old(1);
            end
        end
        
        function prior = PriorPm1(prior)
            % L1 optimization : p(x) ~ 0.5 * delta(x - 1) + 0.5 * delta(x + 1)
            prior.av_mess_old = prior.av_mess;
            prior.av_mess = tanh(prior.R ./ prior.S2);
            prior.var_mess = max(1e-18, 1 - prior.av_mess.^2);
            %             prior.var_mess = 1 - prior.av_mess.^2;
        end
        
        function prior = PriorComplex(prior)
            % Complex case
            var_gauss = prior.param_2; mean_gauss = prior.param_1; R_ = prior.R; S2_ = prior.S2; rho_ = prior.rho;
            prior.av_mess_old = prior.av_mess;
            s2 = var_gauss .* S2_./ (var_gauss + S2_);
            M = (var_gauss .* R_ + S2_ .* mean_gauss) ./ (S2_ + var_gauss);
            alpha_ = abs(mean_gauss).^2 ./ var_gauss + abs(R_).^2 ./ S2_ - abs(M).^2 ./ s2;
            Z = (1 - rho_) .* exp(-0.5 .* abs(R_).^2 ./ S2_) + rho_ .* S2_ ./ (S2_ + var_gauss) .* exp(-alpha_ ./ 2);
            prior.av_mess = (rho_ .* S2_ ./ (S2_ + var_gauss) .* M .* exp(-alpha_ ./ 2) ) ./ Z;
            prior.var_mess = max(1e-18, (rho_ ./ Z .* S2_ ./ (S2_ + var_gauss) .* (2 .* s2 + abs(M).^2) .* exp(-alpha_ ./ 2) - abs(prior.av_mess).^2) ./ 2);
            prior.var_mess(~isfinite(prior.var_mess) ) = 0;            
        end
        
    end
    
end


