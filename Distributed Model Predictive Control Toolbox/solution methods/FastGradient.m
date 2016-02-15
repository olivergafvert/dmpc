classdef FastGradient < DualDecomposition
    % Implements Nesterovs accelerated gradient method [1] to solve the 
    % distributed MPC problem using dual decomposition. The FastGradient
    % class is a subclass of the DualDecomposition class.
    % 
    % The dual variable updates are performed in the coordinators (see help
    % DualDecomposition for an explanation) which are objects of the
    % GDCoordinator class.
    % 
    % To create a FastGradient object the following need to be specified:
    % 
    % FastGradient(sys, maxIter, tol)
    % FastGradient(sys, maxIter, tol, local_qp_solver)
    % FastGradient(sys, maxIter, tol, local_qp_solver, log)
    % FastGradient(sys, maxIter, tol, local_qp_solver, log, method)
    % FastGradient(sys, maxIter, tol, local_qp_solver, log, method, stepsize)
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
    % 'method' - sets the method to use for the dual variable update. The
    % possible methods are Nesterov's accelerated gradient [1]
    % ('fastgradient') and gradient ascent ('gradientascent').
    % 
    % 'stepsize' - sets the step size. If no step size is chosen the step
    % size is chosen as 1/L where L is the Lipschitz constant of the dual
    % gradient.
    %
    % The FastGradient object is the created as follows, given the DmpcSys object, sys:
    % 
    % sol = FastGradient(sys, 1000, 0.001)
    % 
    % which creates an FastGradient solution object with maxIter = 1000, tol =
    % 0.001 and which uses Matlabs quadprog solver to solve the local 
    % quadratic programming problems.
    % 
    % The following line solves the distributed MPC problem for one time
    % step
    % 
    % sol = sol.step()
    % 
    % [1] - Y. Nesterov, "Introductory lectures on convex optimization",
    %       2004, Springer.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        L; % Lipschitz constant of the dual gradient
        method;
    end
    methods
        function obj = FastGradient(sys, maxIter, thresh, local_qp_solver, log, method, stepsize)
            n = nargin;
            if n < 7
                stepsize = 0;
                if n < 6
                    method = 'fastgradient';
                    if n < 5
                        log = Log.NONE;
                        if n < 4
                            local_qp_solver = Solver.QUADPROG;
                        end
                    end
                end
            end
            obj@DualDecomposition(maxIter, thresh, local_qp_solver, log); %Call super class constructor
            obj = obj.mergeSys(sys);
            if stepsize == 0
                obj.L = obj.calcLipschitzConstant();
                stepsize = 1/obj.L;
            end
            obj.method = method;
            for i=1:obj.n_coordinators
                obj.coordinators{i} = obj.coordinators{i}.setStepSize(stepsize);
                obj.coordinators{i} = obj.coordinators{i}.setMethod(obj.method);
            end
            
        end
        
        function obj = setCoordinator(obj, index, group)
            obj.coordinators{index} = GACoordinator(group);
        end
        
        function L = calcLipschitzConstant(obj)
            % Calculates the Lipschitz constant as described in [1]
            % 
            % L = || E * H^-1 * E^T ||,
            % 
            % where E is the dual gradient and H is the Hessian.
            % 
            % [1] - S. Richter, M. Morari and C. N. Jones, "Towards 
            %       Computational Complexity Certification for Constrained
            %       MPC Based on Lagrange Relaxation and the Fast Gradient
            %       Method", 2011, 50th IEEE Conference on Decision and 
            %       Control and European Control Conference (CDC-ECC).
            [EE, H] = obj.getDualGradient();
            E = EE / H;
            try
                L = max(eig(E*EE'));
                L = round(100000000*L)/100000000;
            catch e
                disp(e)
                L = 1;
            end
            disp(['Lipschitz constant: L = ' num2str(L)])
            disp(['Choosing step size: 1/L = ' num2str(1/L)])
        end
        
        function name = getName(obj)
            name = 'FastGradient';
        end
    end
end
