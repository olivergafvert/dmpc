classdef MultiStep < DualDecomposition
    % Implements a multi-step gradient method to solve the 
    % distributed MPC problem using dual decomposition. The MultiStep
    % class is a subclass of the DualDecomposition class.
    % 
    % If we denote the dual variable by lambda the multi-step gradient
    % method performes dual variable updates according to
    % 
    % lambda^{+} = lambda + alpha*(d(lambda)) + beta(lambda - lambda_prev)
    % 
    % where alpha and beta are step sizes (tunable parameters) and
    % d(lambda) is the dual gradient.
    % 
    % The dual variable updates are performed in the coordinators (see help
    % DualDecomposition for an explanation) which are objects of the
    % MSCoordinator class.
    % 
    % To create a MultiStep object the following need to be specified:
    % 
    % MultiStep(sys, maxIter, tol)
    % MultiStep(sys, maxIter, tol, local_qp_solver)
    % MultiStep(sys, maxIter, tol, local_qp_solver, log)
    % 
    % 'sys' - a DmpcSys object describing a distributed system of MPC
    % controllers
    % 
    % 'maxIter' - the maximum number of iterations (i.e updates of the dual
    % variables) that is allowed.
    % 
    % 'tol' - the maximum allowed differance (max(abs(\lambda_1 - \lambda_2))
    % in the dual variables between the subsystems.
    % 
    % 'local_qp_solver' - the quadratic programming problem solver that
    % should be used to solve the local quadratic programming problems
    % (Default is Matlabs quadprog).
    % 
    % 'log' - determines what should be logged using a Log object. 'log'
    % can be either a Log object or one of the macros Log.* (see help Log).
    %
    % The MultiStep object is the created as follows, given the DmpcSys object, sys:
    % 
    % sol = MultiStep(sys, 1000, 0.001)
    % 
    % which creates an MultiStep solution object with maxIter = 1000, tol =
    % 0.001 and which uses Matlabs quadprog solver to solve the local 
    % quadratic programming problems.
    % 
    % The following line solves the distributed MPC problem for one time
    % step
    % 
    % sol = sol.step()
    
    %   Author: Oliver Gäfvert
    %   Email: oliverg@kth.se
    properties
        alpha; % Step size of the multi-step gradient method
        beta; % Step size of the multi-step gradient method
    end
    methods
        function obj = MultiStep(sys, maxIter, thresh, local_qp_solver, log)
            if nargin < 5
                log = Log.NONE;
                if nargin < 4
                    local_qp_solver = Solver.QUADPROG;
                end
            end
            obj@DualDecomposition(maxIter, thresh, local_qp_solver, log); %Call super class constructor
            obj = obj.mergeSys(sys);
            [obj.alpha, obj.beta] = obj.calcStepSizes();
            for i=1:obj.n_coordinators
                obj.coordinators{i} = obj.coordinators{i}.setStepSize(obj.alpha, obj.beta);
            end
        end
        
        function obj = setCoordinator(obj, index, group)
            obj.coordinators{index} = MSCoordinator(group);
        end
        
        function [alpha, beta] = calcStepSizes(obj)
            % Calculates the step sizes as described in [1] for strongly
            % convex problems
            % 
            % alpha = ( 2 /( sqrt(lambda_n(E * H^-1 * E^T)) + 
            %                sqrt(lambda_1(E * H^-1 * E^T)) ) )^2
            % beta = ( ( sqrt(lambda_n(E * H^-1 * E^T)) - 
            %            sqrt(lambda_1(E * H^-1 * E^T) ) /
            %          ( sqrt(lambda_n(E * H^-1 * E^T)) + 
            %            sqrt(lambda_1(E * H^-1 * E^T) ) )^2
            % 
            % where E is the dual gradient, H is the Hessian, lambda_1
            % denotes the smallest non-zero eigenvalue and lambda_n denotes
            % the largest eigenvalue.
            % 
            % [1] - Euhanna Ghadimi, Imam Shames, Mikael Johansson, 
            %       "Multi-step gradient methods for networked 
            %       optimization", 2013.
            [EE, H] = obj.getDualGradient();
            E = EE / H;
            try
                sqmax = sqrt(max(eig(E*EE')));
                sqmin = sqrt(find(eig(E*EE')>0, 1));
            catch e
                disp(e)
                sqmax = 1;
                sqmin = 1;
            end

            alpha = (2 / (sqmax + sqmin))^2;
            beta = ((sqmax-sqmin)/(sqmax+sqmin))^2;

            disp(['Choosing alpha = ' num2str(alpha)])
            disp(['Choosing beta = ' num2str(beta)])
        end
        
        function name = getName(obj)
            name = 'MultiStep';
        end
    end
end
