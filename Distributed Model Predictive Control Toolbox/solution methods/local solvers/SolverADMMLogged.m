classdef SolverADMMLogged < SolverADMM
    % An ADMM solver for a centralized QP.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        logger;
        u_format;
        z_format;
        k; %time step
        iterations; 
    end
    
    methods
        function obj = SolverADMMLogged(H, A, b, Bx, x_k, maxIter, tol, ext)
            % See help SolverADMM
            obj@SolverADMM(H, A, b, Bx, x_k, maxIter, tol);
            obj.logger = Log([ext 'Solver'], Log.SOLVER);
            obj = obj.initLog();
            obj.k = 0;
            obj.iterations = 0;
        end
        
        function obj = update(obj, x_k)
            % Moves the solver to the next time-step.
            obj.c = obj.b + obj.Bx*x_k;
            if ~isempty(obj.z)
                obj.z(1:(end-1), :) = obj.z(2:end, :);
            end
            if ~isempty(obj.u)
                obj.u(1:(end-1), :) = obj.u(2:end, :);
            end
            obj.iterations = 0;
            obj.k = obj.k + 1;
        end
        
        function [obj, x, flag] = solve(obj, q, x)
            % Solves the QP using optimal step size and optimal
            % over-relaxation constant according to [1].
            % 
            % [1] - Euhanna Ghadimi, Andr� Teixeira, Iman Shames, Mikael 
            % Johansson, Optimal parameter selection for the alternating 
            % direction method of multipliers (ADMM): quadratic problems, 
            % "http://arxiv.org/abs/1306.2454", 2013.
            obj.iterations = obj.iterations + 1;
            iter = 0;
            u = obj.u;
            z = obj.z;
            
            AQ = obj.A*obj.Q;
            Qq = AQ*q;
            Qa = AQ*obj.rhoA;
            A_aug = Qq-Qa*obj.c;
            
            
            x = obj.Q * (q + obj.rhoA*(z+u-obj.c));
            Ax = obj.A*x;
            z = max(0, -Ax+obj.c -u);
            u = u + Ax-obj.c + z;
            
            z_p = z;
            Ax = A_aug + Qa*(z+u);
            z = max(0, -obj.alpha*(Ax-obj.c) + (1-obj.alpha)*z-u);
            r = obj.alpha*(Ax+z-obj.c) + (1-obj.alpha)*(z-z_p);
            s = obj.rhoA*(z-z_p);
            u = u + r;
            
            cache_size = 100;
            cache_z = zeros(cache_size, length(z));
            cache_u = zeros(cache_size, length(u));
            cache_norms = zeros(cache_size, 2);
            c_iter = 0;
            s_n=10;r_n=10;
            while iter < obj.maxIter && max(s_n, r_n) > obj.thresh
                iter = iter+1;
                c_iter = c_iter+1;
                z_p = z;
                Ax = A_aug + Qa*(z+u);
                z = max(0, -obj.alpha*(Ax-obj.c) + (1-obj.alpha)*z-u);
                r = obj.alpha*(Ax+z-obj.c) + (1-obj.alpha)*(z-z_p);
                s = obj.rhoA*(z-z_p);
                u = u + r;
                s_n = norm(s, 2); r_n = norm(r, 2);
                
                cache_z(c_iter, :) = z';
                cache_u(c_iter, :) = u';
                cache_norms(c_iter, :) = [s_n r_n];
                if mod(iter, cache_size) == 0
                    obj = obj.logCache(cache_z, cache_u, cache_norms, iter, c_iter, obj.k);
                    c_iter = 0;
                end
                    
            end
            x = obj.Q * (q + obj.rhoA*(z+u-obj.c));
            obj.z = z;
            obj.u = u;
            if iter == obj.maxIter 
                flag = -2;
            else
                flag = 1;
            end
            obj = obj.logCache(cache_z(1:c_iter, :), cache_u(1:c_iter, :), cache_norms(1:c_iter, :), iter, c_iter, obj.k);
            obj.logger = obj.logger.saveToFile();
            %disp(['Inner ADMM iterations: ' num2str(iter)]);
        end
        
        % Logging functions
        function obj = initLog(obj)
            % Log information about this object
            data = [];
            obj.logger = obj.logger.log(data);
            f = '%12.6f';
            zf = ['%12d' f];
            for i=2:size(obj.A, 1)
                zf = [zf f];
            end
            obj.z_format = [zf '\n'];
            obj.u_format = obj.z_format;
        end
        
        function obj = logCache(obj, cache_z, cache_u, cache_norms, iter, offset, time_step)
            z_header = sprintf('%12s%12s\n', 'Iteration', 'z');
            z_str = sprintf(obj.z_format, [(iter-offset+1):iter;cache_z']);
            u_header = sprintf('%12s%12s\n', 'Iteration', 'u');
            u_str = sprintf(obj.u_format, [(iter-offset+1):iter;cache_u']);
            norms_header = sprintf('%12s %12s %12s\n', 'Iteration', 's-norm', 'r-norm');
            norms = sprintf('%12d %12.6f %12.6f\n', [(iter-offset+1):iter;cache_norms']);
            data = ['\n--------------------------------------------\n Solver Iteration: ' ...
                num2str(obj.iterations) '\n Time step k = ' num2str(time_step) '\n\n z: \n' ...
                z_header z_str '\n u:\n' u_header u_str '\n' norms_header norms '\n'];
            obj.logger = obj.logger.log(data);
        end
        
    end
end
            