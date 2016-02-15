classdef DualDecomposition
    % The DualDecomposition class is a solution object in the sense that
    % given a DmpcSys object it solves the global MPC problem described by
    % the DmpcSys object in a distributed way using dual decomposition. The
    % DualDecomposition class is meant to serve as a super class for 
    % methods that use dual decomposition to solve the distributed MPC
    % problem. It is in itself not a complete implementation of dual
    % decomposition since it has no implemented dual variable update.
    % 
    % The DualDecomposition class (or subclasses) uses the DualObject
    % class (or subclasses) to represent subsystems and the Coordinator
    % class (or subclasses) to manage the communication. To simulate 
    % communication in the distributed system (i.e message passing between 
    % objects and coordinators) an object of the class Communication is used.
    % 
    % To create create a DualDecomposition object the following parameters
    % need to be specified:
    % 
    % DualDecomposition(maxIer, tol)
    % DualDecomposition(maxIer, tol, local_qp_solver)
    % DualDecomposition(maxIer, tol, local_qp_solver, log)
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
    % The DualDecomposition object is then cerated as follows:
    % 
    % sol = DualDecomposition(1000, 0.001)
    % 
    % This line creates a DualDecomposition object with maxIter = 1000, tol
    % = 0.001 and that uses Matlabs quadprog solver for the local quadtratic
    % programming problems. The following line maps a DmpcSys object, sys, to
    % this solution object:
    % 
    % sol = sol.mergeSys(sys)
    % 
    % now the following solves the distributed MPC problem for one time step
    % 
    % sol = sol.step()
    %
    % and the following simulates 10 time steps
    % 
    % sol = sol.sim(10)
    
    %   Author: Oliver Gäfvert
    %   Email: oliverg@kth.se
    properties
        com; % The Communication object that simulates comminucation between subsystems
        histStates; % A matrix of previous states
        histInputs; % A matrix of previous inpus signals
        histIter; % A vector of the previous number of iterations until convergence
        
        objects; % A cell array of the objects (subsystems) of this solution object
        n_objects; % the number of objects (subsystems)
        object_map; % maps object ids to positions in the 'objects' cell array
        
        coordinators; % A cell array of the coordinators (groups) of this solution object
        n_coordinators; % the number of coordinators (groups)
        coordinator_map; % maps coordinator ids to positions in the 'coordinators' cell array
        
        maxIter; % The maximum number of iterations (dual variable updates) allowed
        tol; % the maximum allowed L_2 norm differance in the dual variables between the subsystems.
        
        local_qp_solver; % the quadratic programming problem solver that
                         % should be used to solve the local quadratic programming problems
                         % (Default is Matlabs quadprog).
        logger; % a Log object that creates a log file.
    end
    
    methods
        function obj = DualDecomposition(maxIter, tol, local_qp_solver, log)
            % See help DualDecomposition
            obj.com = Communication();
            obj.objects = cell(1, 1);
            obj.n_objects = 0;
            obj.coordinators = cell(1, 1);
            obj.n_coordinators = 0;
            obj.histIter = [];
            obj.maxIter = maxIter;
            obj.tol = tol;
            if nargin < 4
                log = Log.NONE;
                if nargin < 3
                    local_qp_solver = 'quadprog';
                end
            end
            if strcmp(log, Log.NONE)
                obj.logger = 0;
            else
                if isa(log, 'Log')
                    obj.logger = log;
                else
                    obj = obj.initLog(log);
                end
            end
            obj.local_qp_solver = local_qp_solver;
        end
        
        function obj = mergeSys(obj, sys)
            % obj = obj.merge(sys) takes
            % sys - a DmpcSys object 
            % as input and merges the DmpcSubsystem objects into
            % DualObject (or subclasses) objects and the DmpcGroup objects
            % into Coordinator (or subclasses) objects.
            obj.object_map = zeros(2, 1);
            for i=1:sys.n_subsystems
                obj.n_objects = obj.n_objects + 1;
                obj = obj.setObject(obj.n_objects, sys.subsystems{i});
                obj.object_map(obj.objects{obj.n_objects}.id) = obj.n_objects;
            end
            for i=1:sys.n_groups
                obj.n_coordinators = obj.n_coordinators + 1;
                obj = obj.setCoordinator(obj.n_coordinators, sys.groups{i});
                obj.coordinator_map(obj.coordinators{obj.n_coordinators}.id) = obj.n_coordinators;
            end
            
            if obj.shouldLog()
                for i=1:sys.n_subsystems
                    obj.objects{i} = obj.objects{i}.initLog(obj.logger.type, obj.getName());
                end
                for i=1:sys.n_groups
                    obj.coordinators{i} = obj.coordinators{i}.initLog(obj.logger.type, obj.getName());
                end
            end
        end
        
        function obj = setCoordinator(obj, index, group)
            obj.coordinators{index} = Coordinator(group);
        end
        
        function obj = setObject(obj, index, object)
            obj.objects{index} = DualObject(object);
            obj.objects{index} = obj.objects{index}.initSolver(obj.local_qp_solver);
        end
        
        function pos = getCoordinatorPosition(obj, id)
            pos = obj.coordinator_map(id);
        end
        
        function pos = getObjectPosition(obj, id)
            pos = obj.object_map(id);
        end
        
        function obj = step(obj, update, updateSolver)
            % Solves the distributed MPC problem for one time step.
            if nargin < 2
                update = 1;
                if nargin < 3
                    updateSolver = 1;
                end
            end
            epsilonCoord = ones(obj.n_coordinators, 1);
            thCoord = ones(obj.n_coordinators, 1)*obj.tol;
            epsilonObj = ones(obj.n_objects, 1);
            thObj = ones(obj.n_objects, 1)*obj.tol;
            iter = 0;
            ids_obj = 1:obj.n_objects;
            ids_coordinators = 1:obj.n_coordinators;
            
            obj.com = obj.com.clearOut();
            while (sum(epsilonCoord > thCoord) > 0 || sum(epsilonObj > thObj) > 0) && iter < obj.maxIter
                iter = iter + 1;
                c_objects = obj.objects;
                out = cell(obj.n_objects, 1);
                
                for i=ids_obj
                    out{i} = obj.com.getMessageDMPC(obj.objects{i}.id);
                end
                
                temp_message = cell(obj.n_objects, 1);

                parfor i=ids_obj
                    [c_objects{i}, epsilonObj(i)] = c_objects{i}.setDualVariable(out{i});
                    [c_objects{i}, message] = c_objects{i}.step();
                    temp_message{i} = message;
                end
                obj.objects = c_objects;
                for i=ids_obj
                    obj.com = obj.com.toGroup(obj.objects{i}.id, temp_message{i});
                end
                
                obj.com = obj.com.clearOut();
                
                for i=ids_coordinators
                    [obj.coordinators{i}, message, epsilonCoord(i)] = ...
                        obj.coordinators{i}.evalCoupledVariables(obj.com.getMessageGroup(obj.coordinators{i}.id));
                    obj.com = obj.com.toDMPC(obj.coordinators{i}.id, message);
                end
                obj.com = obj.com.clearIn();
            end
            
            disp(['Iterations: ' num2str(iter)])
            
            obj.histIter = [obj.histIter, iter];
            
            if update == 1
                for i=1:obj.n_coordinators
                    obj.coordinators{i} = obj.coordinators{i}.update();
                end
            end
            
            for i=1:obj.n_objects
                obj.objects{i} = obj.objects{i}.setDualVariable(obj.com.getMessageDMPC(obj.objects{i}.id));
                if update == 1
                    obj.objects{i} = obj.objects{i}.updateState();
                end
                if updateSolver == 1
                    obj.objects{i} = obj.objects{i}.updateSolver();
                end
            end
            
            if obj.shouldLog()
                obj = obj.flushLog();
            end
        end
        
        function obj = initLog(obj, log)
            obj.logger = Log(obj.getName(), log);
            
        end
        
        function b = shouldLog(obj)
            if isa(obj.logger, 'Log')
                b = 1;
                return;
            end
            b=0;
        end
        
        function obj = flushLog(obj)
            for i=1:obj.n_objects
                obj.objects{i}.logger = obj.objects{i}.logger.saveToFile();
            end
            for i=1:obj.n_coordinators
                obj.coordinators{i}.logger = obj.coordinators{i}.logger.saveToFile();
            end
        end
        
        function obj = log(obj)
            
        end
            
        
        function [obj, histIter] = sim(obj, steps, hist)
            % Solve the distributed MPC problem for 'steps' number of time
            % steps using the local LTI objects to simulate future states.
            if nargin < 3
                hist = '';
            end
            obj = obj.startHist();
            for i=1:steps
                obj = obj.step();
                obj = obj.addHist();
            end
            histIter = obj.histIter;
            
            if nargin > 2 && strcmp(hist, 'plot')
                obj.plotHistory();
            end
        end
        
        function [EE, Hessian] = getDualGradient(obj)
            % Returns the dual gradient (and the hessian) of the global 
            % MPC problem.
            Hessian = [];
            ctr = 0;
            map = zeros(10, 2);
            
            for i=1:obj.n_objects
                H = obj.objects{i}.mpc.H;
                map(obj.objects{i}.id, :) = [ctr, (ctr+size(H, 1))];
                ctr = ctr + size(H, 1);
                Hessian = blkdiag(Hessian, H);
            end
            
            EE = [];
            for i=1:obj.n_coordinators
                gid = obj.coordinators{i}.id;
                members = obj.coordinators{i}.members';
                pos = obj.object_map(members(1), 1);
                E1 = obj.objects{pos}.getCouplingMatrix(gid);
                pos = obj.object_map(members(2), 1);
                E2 = obj.objects{pos}.getCouplingMatrix(gid);
                
                EEt = zeros(size(E1, 1), size(Hessian, 2));
                
                pos1 = map(members(1), :);
                pos2 = map(members(2), :);
                EEt(:, (pos1(1)+1):pos1(2)) = E1;
                
                EEt(:, (pos2(1)+1):pos2(2)) = -E2;
                EE = [EE; EEt];
            end
        end
        
        function plotHistory(obj)
            figure
            plot(obj.histIter)
            title('Number of iterations')
        end
        
        function obj = startHist(obj)
            obj.histStates = cell(obj.n_objects, 1);
            obj.histInputs = cell(obj.n_objects, 1);
            obj = obj.addHist();
        end
        
        function obj = addHist(obj)
            for i=1:obj.n_objects
                obj.histStates{i} = [obj.histStates{i}, obj.objects{i}.mpc.x_k];
                obj.histInputs{i} = [obj.histInputs{i}, obj.objects{i}.getPrev_u()];
            end
        end
        
        function name = getName(obj)
            name = 'DualDecomposition';
        end
        
    end
end
