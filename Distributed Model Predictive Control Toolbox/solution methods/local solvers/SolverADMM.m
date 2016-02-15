classdef SolverADMM < Solver
    % An ADMM solver for a centralized QP.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        Q;
        rhoA;
        rho;
        alpha;
        
        z;
        u;
        
        maxIter;
        thresh;
    end
    
    methods
        function obj = SolverADMM(H, A, b, Bx, x_k, maxIter, tol)
            % See help SolverADMM
            obj@Solver(H, [], A, b, Bx, x_k)
            n = nargin;
            if n<7
                obj.thresh = 0.000001;
                if n<6
                    obj.maxIter = 1000000;
                else
                    obj.maxIter = maxIter;
                end
            else
                obj.thresh = tol;
                obj.maxIter = maxIter;
            end
            [obj.rho, obj.alpha] = obj.computeStepSize(obj.H ./ 2, obj.A);
            obj.Q = -(inv(obj.H + obj.rho*(obj.A' * obj.A)));
            obj.rhoA = obj.rho*obj.A';
            obj.z = [];
            obj.u = [];
        end
        
        function [rho, alpha] = computeStepSize(obj, H, A)
            % Computes the optimal step size as in [1]
            % 
            % rho = 1/sqrt( lambda_1(A * H^-1 * A^T) * ...
            %               lambda_n(A * H^-1 * A^T) )
            % 
            % where lambda_1 denotes the smallest (non-zero) eigenvalue 
            % and lambda_n denotes the largest eigenvalue. And sets the
            % optimal over-relaxation parameter, alpha, according to [1].
            % 
            % [1] - Euhanna Ghadimi, Andr� Teixeira, Iman Shames, Mikael 
            % Johansson, Optimal parameter selection for the alternating 
            % direction method of multipliers (ADMM): quadratic problems, 
            % "http://arxiv.org/abs/1306.2454", 2013.
            eigA = eig((A / H) * A');
            minA = find(eigA>0, 1);
            maxA = max(eigA);
            
            rho = 1/sqrt(maxA*minA);
            alpha = 2; % over-relaxation parameter
        end
        
        function obj = update(obj, x_k)
            % Moves the solver to the next time-step.
            obj.c = obj.b + obj.Bx*x_k;
        end
        
        function [obj, x, flag] = solve(obj, q, x)
            % Solves the QP using optimal step size and optimal
            % over-relaxation constant according to [1].
            % 
            % [1] - Euhanna Ghadimi, Andr� Teixeira, Iman Shames, Mikael 
            % Johansson, Optimal parameter selection for the alternating 
            % direction method of multipliers (ADMM): quadratic problems, 
            % "http://arxiv.org/abs/1306.2454", 2013.
            iter = 0;
            u = obj.u;
            z = obj.z;
            
            if isempty(u)
                u=0;
            end
            if isempty(z)
                z = 0;
            end
            
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
            
            while iter < obj.maxIter && max(norm(s, 2), norm(r, 2)) > obj.thresh
                iter = iter+1;
                z_p = z;
                Ax = A_aug + Qa*(z+u);
                z = max(0, -obj.alpha*(Ax-obj.c) + (1-obj.alpha)*z-u);
                r = obj.alpha*(Ax+z-obj.c) + (1-obj.alpha)*(z-z_p);
                s = obj.rhoA*(z-z_p);
                u = u + r;
            end
            x = obj.Q * (q + obj.rhoA*(z+u-obj.c));
            obj.z = z;
            obj.u = u;
            if iter == obj.maxIter 
                flag = -2;
            else
                flag = 1;
            end
           %disp(['Inner ADMM iterations: ' num2str(iter)]);
        end
    end
end
            