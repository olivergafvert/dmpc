classdef ADMMObject < DualObject
    % The ADMMObject class is a dual object within the ADMM solution object
    % (see help DualObject, ADMM). 
    % 
    % ADMMObject(subsys, rho, alpha, solver) creates an ADMMObject using
    % 'subsys' - a DmcpSubsystem object - with step size 'rho',
    % over-relaxation constant 'alpha' and which uses 'solver' to solver
    % the local quadratic programming problems.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        %Dual variables
        z; % the z-variable in ADMM
        rho; % step size of ADMM
        
        %QP arguments
        H;
        
        alpha; % over-relaxation constant
        
        %Coupling matrices
        Ec;
        Mc;
        
        % Logging properties
        z_format;
    end
    
    methods
        function obj = ADMMObject(subsys, rho, alpha, solver, maxIter, tol)
            % See help ADMMObject
            n = nargin;
            if n < 5
                maxIter = 1000000;
                if n < 6
                    tol = 0.001;
                end
            end
            obj@DualObject(subsys);
            obj.iter = 0;
            obj.prev_u = [];
            obj.rho = rho;
            obj.alpha = alpha;
            obj.z = zeros(size(obj.mpc.LTI.c, 1), obj.mpc.Nc);
            obj.dual_variable = zeros(size(obj.mpc.LTI.c, 1)*obj.mpc.Nc, 1);
            [obj.Ec, obj.Mc] = obj.getCouplingMatrixAll();
            obj.H = obj.mpc.H +(obj.rho)*(obj.Ec'*obj.Ec);
            obj.ft =  (obj.rho*(obj.Ec' * obj.Mc * obj.mpc.x_k)) + obj.mpc.F*obj.mpc.x_k;
            obj.f = obj.ft + (obj.Ec' *((obj.dual_variable) - (obj.rho) * obj.z(:)));
            obj = obj.initSolver(solver, obj.H, obj.mpc.F, obj.mpc.Ac, ...
                obj.mpc.b, obj.mpc.Bx, obj.mpc.x_k, maxIter, tol);
        end
        
        function [E, M] = getCouplingMatrixAll(obj)
            % Returns the matrix Ec and Mc so that 
            % E_i*y_i = Ec*u_i + Mc*x_i(k|k)
            rows = zeros(size(obj.dual_variable));
            temp = sign(obj.dual_coord);
            tr = size(temp, 1);
            for i=1:obj.mpc.Nc
                rows(((i-1)*tr+1):(i*tr), 1) = temp;
            end
            E = zeros(size(obj.Conv));
            index = find(rows>0);
            E(index, :) = obj.Conv(index, :);
            M(index, :) = obj.M(index, :);
            if isempty(E)
                E=1;
            end
            if isempty(M)
                M=1;
            end
        end
        
        function y = getConstY(obj)
            % Returns E_i * y_i (i.e a vector that is coupled to
            % other systems by E_i*y_i = z_i)
            y = obj.Ec*obj.prev_u + obj.Mc(1:size(obj.Ec, 1), :)*obj.mpc.x_k;
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
            if ~isempty(obj.z)
                tmp_dual = reshape(obj.dual_variable, [size(obj.mpc.LTI.c, 1) obj.mpc.Nc]);
                for i=1:size(obj.mpc.LTI.c, 1)
                    obj.z(i, 1:(end-1)) = obj.z(i, 2:end);
                    obj.z(i, end) = 2*obj.z(i, end-1) - obj.z(i, end-2);
                    tmp_dual(i, 1:(end-1)) = tmp_dual(i, 2:end);
                    tmp_dual(i, end) = 2*tmp_dual(i, end-1) - tmp_dual(i, end-2);
                end
                obj.dual_variable = tmp_dual(:);
            end
            obj.iter = 0;
            obj.k = obj.k + 1;
            obj.ft = ((obj.rho)*(obj.mpc.x_k' * obj.Mc' *obj.Ec))' + obj.mpc.F*obj.mpc.x_k;
            obj.f = obj.ft + (obj.Ec' *((obj.dual_variable) - (obj.rho) * obj.z(:)));
        end
        
        function obj = updateSolver(obj)
            obj.ft = ((obj.rho)*(obj.mpc.x_k' * obj.Mc' *obj.Ec))' + obj.mpc.F*obj.mpc.x_k;
            obj.solver = obj.solver.update(obj.mpc.x_k);
        end
        
        function message = sendMessage(obj)
            % Prepares a message to send to Coordinator object(s). The
            % messages to each Coordinator contains the value 
            % E_{ij}*y + lambda/rho where lambda is the dual variable and
            % rho is the step size.
            index_vector = find(obj.dual_coord > 0);
            ys = obj.getY();
            ys = obj.alpha*ys + (1-obj.alpha)*obj.z(:) + obj.dual_variable ./ obj.rho;
            ys = reshape(ys, [size(obj.z, 1), obj.mpc.Nc]);
            message = [obj.dual_coord(index_vector) , ys(index_vector, :)];
        end
        
        function [obj, epsilon] = setDualVariable(obj, message)
            % Sets the dual variable (lambda) for this object (subsystem).
            % The message should contain the updated 'z' varialbe and
            % using this the updated dual variable is computed as
            % 
            % lambda = lambda + rho*(E*y - z)
            % 
            rows = size(message, 1);
            if rows == 0
                epsilon = 100;
                return 
            end
            z_temp = obj.z(:);
            for i=1:rows
                obj.z(obj.coord_dual(message(i, 1)), :) = message(i, 2:(end));
            end
            z = obj.z(:);
            y = obj.alpha*obj.getConstY() + (1-obj.alpha)*z_temp;
            grad = y - z;
            ep = norm(obj.rho*(z-z_temp), 2);
            epsilon = norm(grad, 2);
            epsilon = max(epsilon, ep);
            %disp(['Epsilon: ' num2str(epsilon)])
            
            obj.dual_variable = obj.dual_variable + obj.rho*(grad);
            obj.f = obj.ft + (obj.Ec' *((obj.dual_variable) - (obj.rho) * obj.z(:)));
        end
        
        % Logging functions
        function obj = initLog(obj, type, ext)
            if nargin < 3
                ext = '';
            end
            if strcmp(type, Log.ALL) || strcmp(type, Log.COORDINATORS)
                obj.logger = Log([ext 'Subsystem' num2str(obj.id)], type);
                c_title = sprintf('%12s %12s\n', 'Coord id', 'y-index');
                index_vector = find(obj.dual_coord > 0);
                coupling_str = sprintf('%12d %12d\n', [obj.dual_coord(index_vector)  index_vector]');
                % Log information about this object
                data = ['Log of Subsystem' num2str(obj.id) '\n'...
                ' x_0 : ' mat2str(obj.mpc.x_k') '\n\n Couplings:\n' c_title '\n' coupling_str '\n'];
                obj.logger = obj.logger.log(data);
                f = '%12.6f';
                uf = f;
                for i=2:size(obj.mpc.LTI.b, 2)
                    uf = [uf f];
                end
                obj.u_format = [uf '\n'];
                zf = f;
                for i=2:size(obj.z, 1)
                    zf = [zf f];
                end
                obj.z_format = [zf '\n'];
                obj.l_format = obj.z_format;
                obj.y_format = obj.l_format;
            end
        end
        
        function obj = log(obj)
            u = reshape(obj.prev_u(1:(size(obj.mpc.LTI.b, 2)*obj.mpc.Nc)), [size(obj.mpc.LTI.b, 2), obj.mpc.Nc]);
            u_str = sprintf(obj.u_format, u);
            l_str = sprintf(obj.l_format, reshape(obj.dual_variable, [size(obj.z, 1), obj.mpc.Nc]));
            z_str = sprintf(obj.z_format, obj.z);
            ys = obj.getY();
            ys = reshape(ys, [size(obj.z, 1), obj.mpc.Nc]);
            y_str = sprintf(obj.y_format, ys);
            data = ['\n--------------------------------------------\nIteration: '...
                num2str(obj.iter) ' \nTime step k = ' num2str(obj.k) '\n\n y (t\\index) : \n' ...
                y_str '\n u (t\\index) :\n'...
                u_str '\n x : '...
                mat2str(obj.mpc.x_k') '\n\n dual variable (t\\index) : \n' ...
                l_str '\n\n z (t\\index) :\n' z_str '\n'];
            obj.logger = obj.logger.log(data);
        end
    end
end