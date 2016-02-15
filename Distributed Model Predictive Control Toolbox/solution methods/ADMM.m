classdef ADMM < DualDecomposition
    % Implements the Alternating Direction Method of Multipliers (ADMM)
    % using edge variables to solve the distributed MPC problem. The ADMM
    % class is a subclass of the DualDecomposition class.
    % 
    % To create a ADMM object the following need to be specified:
    % 
    % ADMM(sys, maxIter, tol)
    % ADMM(sys, maxIter, tol, local_qp_solver)
    % ADMM(sys, maxIter, tol, local_qp_solver, rho, alpha)
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
    % 'rho' - step size of the ADMM (Default = 1)
    % 
    % 'alpha' - over-relaxation parameter (Default = 1.6)
    %
    % The ADMM object is the created as follows, given the DmpcSys object, sys:
    % 
    % sol = ADMM(sys, 1000, 0.001)
    % 
    % which creates an ADMM solution object with maxIter = 1000, tol =
    % 0.001 and which uses Matlabs quadprog solver to solve the local 
    % quadratic programming problems.
    % 
    % The following line solves the distributed MPC problem for one time
    % step
    % 
    % sol = sol.step()
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        rho; % The step size of the ADMM algorithm
        alpha; %over-relaxation paramter
    end
    methods
        function obj = ADMM(sys, maxIter, tol, local_qp_solver, log, rho, alpha)
            % See help ADMM
            n = nargin;
            if n < 7
                alpha = 1.33;
                if n < 6
                    rho = 1;
                    if n < 5
                        log = Log.NONE;
                        if n < 4
                            local_qp_solver = Solver.QUADPROG;
                        end
                    end
                end
            end
            obj@DualDecomposition(maxIter, tol, local_qp_solver, log); %Call super class constructor
            obj.rho = rho;
            obj.alpha = alpha;
            obj = obj.mergeSys(sys);
        end
        
        function obj = setCoordinator(obj, index, group)
            obj.coordinators{index} = ADMMCoordinator(group);
        end
        
        function obj = setObject(obj, index, object)
            obj.objects{index} = ADMMObject(object, obj.rho, obj.alpha, obj.local_qp_solver);
        end
        
        function name = getName(obj)
            name = 'ADMM';
        end
    end
end
