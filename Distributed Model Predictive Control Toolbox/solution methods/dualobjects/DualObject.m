classdef DualObject
    % A DualObject is a representation of a DmpcSubsystem object within 
    % the solution object that is used to solve the distributed MPC
    % problem. 
    % 
    % DualObject(subsys) creates a DualObject from 'subsys' - a
    % DmpcSubsystem object.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        id; % the id of this object within the distributed system
        mpc; %the Mpc object belonging to this object
        k; %time-step variable
        iter; % the current number of iterations (i.e the number of dual 
              % variable updates that has been performed)
        coordinators; %pointers to other dmpc objects that are coupled to this subsystem
        n_coordinators; % the number of groups coupled with the object
        solver; % the solver for the local quadratic programming problem
        prev_u; % the current solution to the local qp-problem
        dual_variable;
        coord_dual; % a map mapping coordinator to dual variable
        dual_coord; % a map mapping dual variable to coordinator
        
        % The 'Conv' and 'M' matrices are the matrices such that 
        % y = Conv*prev_u + M*mpc.x_k where y is the output vector of the
        % local system.
        Conv; 
        M;
        
        % QP-argument 
        f;
        ft;
        
        % Logging properties
        logger;
        u_format;
        l_format;
        y_format;
    end
    
    methods
        function obj = DualObject(subsys)
            % See help DmpcSubsystem
            obj.id = subsys.id;
            obj.k = 0;
            obj.iter = 0;
            obj.mpc = subsys.mpc;
            obj.n_coordinators = 0;
            obj.coordinators = zeros(1);
            
            obj.dual_variable = zeros(size(obj.mpc.LTI.c, 1), obj.mpc.Nc);
            obj.dual_coord = zeros(size(obj.dual_variable, 1), 1);
            obj.coord_dual = zeros(1);
            
            if subsys.n_groups ~= 0
                for i=1:subsys.n_groups
                    obj = obj.addCoordinator(subsys.groups(i), ...
                        subsys.getCoupledVar(subsys.groups(i)));
                end
            end
            
            obj = obj.initConvMatrix();
            obj.logger = 0;
            obj.ft = obj.mpc.F*obj.mpc.x_k;
            obj.f = obj.ft + (obj.Conv' * obj.dual_variable(:));
        end
        
        function obj = addCoordinator(obj, group, y_index)
            % See DmpcSys.addGroup
            obj.n_coordinators = obj.n_coordinators + 1;
            obj.coordinators(obj.n_coordinators) = group;
            obj.coord_dual(group) = y_index;
            obj.dual_coord(y_index) = group;
        end
        
        function obj = initSolver(obj, solver, H, f, A, b, Bx, x, maxIter, tol)
            % Initializes the solver that should be used for the local
            % qp-problems. (Default is Matlabs quadprog function)
            if nargin < 3
                H = obj.mpc.H;
                f = obj.mpc.F*obj.mpc.x_k;
                A = obj.mpc.Ac;
                b = obj.mpc.b;
                Bx = obj.mpc.Bx;
                x = obj.mpc.x_k;
            end
            if nargin < 9
                maxIter = 1000000;
                tol = 0.000001;
            end
            if strcmp(solver, Solver.ADMM)
                obj.solver = SolverADMM(H, A, b, Bx, x, maxIter, tol);
            elseif strcmp(solver, Solver.ADMM_LOGGED)
                 obj.solver = SolverADMMLogged(H, A, b, Bx, x, maxIter, tol, obj.getName());
            else
                obj.solver = Solver(H, f*x, A, b, Bx, x);
            end
        end
        
        function E = getCouplingMatrix(obj, gid)
            % Returns the matrix E encoding the coupling E_{i} * u = z_{i} 
            % that is represented by the group with group id 'gid'.
            yc = find(obj.dual_coord == gid);
            s = size(obj.mpc.LTI.c, 1);
            rows = yc:s:(s*obj.mpc.Nc);
            E = obj.Conv(rows, :);
        end
        
        function obj = initConvMatrix(obj)
            % Initializes the matrices obj.Conv and obj.M that are the 
            % matrices such that 
            % y = obj.Conv*obj.prev_u + obj.M*obj.mpc.x_k, where y is the 
            % output vector of the local system.
            C = obj.mpc.LTI.c;
            CC = C;
            D = obj.mpc.LTI.d;
            DD = D;
            
            for i=1:(obj.mpc.Nc-1)
                CC = blkdiag(CC, C);
                DD = blkdiag(DD, D);
            end
            obj.Conv = CC*obj.mpc.Conv + DD;
            obj.M = CC*obj.mpc.M;
            if obj.mpc.soft > 0
                n = size(obj.mpc.H, 1) - size(obj.Conv, 2);
                obj.Conv = [obj.Conv, zeros(size(obj.Conv, 1), n)];
                obj.M = [obj.M; zeros(n, size(obj.M, 2))];
            end
        end
        
        function [obj, message] = step(obj)
            % Solves the local qp-problem for one time-step.
            obj.iter = obj.iter + 1;
            
            %Solve qp
            [obj.solver, obj.prev_u, flag] = obj.solver.solve(obj.f, obj.prev_u);

            if flag ~= 1 % Check if a solution is found
                disp('Found no solution, returning approximation')
            end
            
            if obj.shouldLog()
                obj = obj.log();
            end
            
            message = obj.sendMessage();
        end
        
        function message = sendMessage(obj)
            % Prepares a message to send to Coordinator object(s). The
            % messages to each Coordinator contains the value 
            % subset of the output vector of this subsystem that is coupled
            % to the Coordinator.
            index_vector = find(obj.dual_coord > 0);
            ys = obj.getY();
            ys = reshape(ys, [size(obj.dual_variable, 1), obj.mpc.Nc]);
            message = [obj.dual_coord(index_vector) , ys(index_vector, :)];
        end
        
        function [obj, epsilon] = setDualVariable(obj, message)
            % Sets the dual variable (lambda) for this object (subsystem).
            % The message should contain subsets of the the updated dual
            % varialbe from each Coordinator, i.e a message should look like
            % 
            % [id of Coordinator j, lambda_ij ;
            %  id of Coordinator k, lambda_ik ;
            %  ...                             ]
            rows = size(message, 1);
            epsilon = 0;
            if rows == 0
                return;
            end
            if isa(message(:, 1), 'cell')
                tmp = zeros(size(rows))
                for i=1:rows
                    m_temp = message(i, 1);
                    tmp(i) = m_temp{1};
                end
                ids = tmp;
            else
                ids = message(:, 1);
            end
            
            for i=1:rows
                obj.dual_variable(obj.coord_dual(ids(i)), :) = message(i, 2:(end));
            end
            obj.f = obj.ft + (obj.Conv' * obj.dual_variable(:));
        end
        
        function y = getY(obj)
            y = obj.Conv*obj.prev_u + obj.M(1:size(obj.Conv, 1), :)*obj.mpc.x_k;
        end
        
        function u = getPrev_u(obj)
            bc = size(obj.mpc.LTI.b, 2);
            if size(obj.prev_u, 1) ~= 0
                u = obj.prev_u(1:bc);
            else
                u = [];
            end
        end
        
        function obj = updateState(obj, x_k)
            % Moves this subsystem to the next time step using the comupted
            % input signals and the LTI model of the system.
            bc = size(obj.mpc.LTI.b, 2);
            if nargin < 2
                obj.mpc.x_k = obj.mpc.LTI.a*obj.mpc.x_k + obj.mpc.LTI.b*obj.prev_u(1:bc);
            else
                obj.mpc.x_k = x_k;
            end
            if ~isempty(obj.prev_u)
                u = obj.prev_u(1:(size(obj.mpc.LTI.b, 2)*obj.mpc.Nc));
                u = reshape(u, [size(obj.mpc.LTI.b, 2) obj.mpc.Nc]);
                for i=1:size(obj.mpc.LTI.b, 2)
                    u(i, 1:(end-1)) = u(i, 2:end);
                    u(i, end) = 2*u(i, end-1) - u(i, end-2);
                end
                u = u(:);
                obj.prev_u(1:(length(u))) = u;
            end
            if ~isempty(obj.dual_variable)
                for i=1:size(obj.mpc.LTI.c, 1)
                    obj.dual_variable(i, 1:(end-1)) = obj.dual_variable(i, 2:end);
                    obj.dual_variable(i, end) = 2*obj.dual_variable(i, end-1) - obj.dual_variable(i, end-2);
                end
            end
            obj.iter = 0;
            obj.k = obj.k + 1;
            obj.ft = obj.mpc.F*obj.mpc.x_k;
            obj.f = obj.ft + (obj.Conv' * obj.dual_variable(:));
        end
        
        function obj = updateSolver(obj)
            obj.solver = obj.solver.update(obj.mpc.x_k);
        end
        
        function print(obj)
            % Prints some information about this object
            disp(['Subsystem id: ' num2str(obj.id)]);
            disp('LTI object information:');
            disp('LTI.a');
            disp(obj.mpc.LTI.a);
            disp('LTI.b');
            disp(obj.mpc.LTI.b);
            disp('LTI.c');
            disp(obj.mpc.LTI.c);
            disp('LTI.d');
            disp(obj.mpc.LTI.d);
        end
        
        function obj = initLog(obj, type, ext)
            if nargin < 3
                ext = '';
            end
            if strcmp(type, Log.ALL) || strcmp(type, Log.SUBSYS)
                obj.logger = Log([ext 'Subsystem' num2str(obj.id)], type);
                c_title = sprintf('%12s %12s\n', 'Coord id', 'y-index');
                index_vector = find(obj.dual_coord > 0);
                coupling_str = sprintf('%12d %12d\n', [obj.dual_coord(index_vector)  index_vector]');
                % Log information about this object
                data = ['Log of Subsystem' num2str(obj.id) '\n'...
                ' x_0 : ' mat2str(obj.mpc.x_k') '\n\n Couplings:\n' coupling_str ];
                obj.logger = obj.logger.log(data);
                f = '%12.6f';
                uf = f;
                for i=2:size(obj.mpc.LTI.b, 2)
                    uf = [uf f];
                end
                obj.u_format = [uf '\n'];
                lf = f;
                for i=2:size(obj.dual_variable, 1)
                    lf = [lf f];
                end
                obj.l_format = [lf '\n'];
                obj.y_format = obj.l_format;
            end
        end
        
        function b = shouldLog(obj)
            if isa(obj.logger, 'Log')
                b = 1;
                return;
            end
            b=0;
        end
        
        function obj = log(obj)
            u = reshape(obj.prev_u(1:(size(obj.mpc.LTI.b, 2)*obj.mpc.Nc)), [size(obj.mpc.LTI.b, 2), obj.mpc.Nc]);
            u_str = sprintf(obj.u_format, u);
            l_str = sprintf(obj.l_format, obj.dual_variable);
            ys = obj.getY();
            ys = reshape(ys, [size(obj.dual_variable, 1), obj.mpc.Nc]);
            y_str = sprintf(obj.y_format, ys);
            data = ['\n--------------------------------------------\nIteration: '...
                num2str(obj.iter) ' \nTime step k = ' num2str(obj.k) '\n\n y (t\\index) : \n' ...
                y_str '\n u (t\\index) :\n'...
                u_str '\n x : '...
                mat2str(obj.mpc.x_k') '\n\n dual variable (t\\index) : \n' ...
                l_str '\n'];
            obj.logger = obj.logger.log(data);
        end
        
        function name = getName(obj)
            name = ['Subsystem' num2str(obj.id)];
        end
        
    end
end