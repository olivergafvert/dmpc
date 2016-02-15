classdef ADMMCoordinator < Coordinator
    % A ADMMCoordinator manages the communication of coupled varaibles and dual
    % variables between two ADMMObjects. A ADMMCoordinator object manages
    % the update of the 'z' variable in the ADMM algorithm.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        ctr;
        z;
    end
    
    methods
        function obj = ADMMCoordinator(group)
            % See help ADMMCoordinator
            obj@Coordinator(group);
            obj.ctr = 0;
        end
        
        function obj = updateCoupledVariables(obj, id, val)
            % Stores the (y + lambda/rho) value for each of this Coordinators
            % members. When all the members has sent a new value of their
            % of this variable those values are used to compute the 'z'
            % varaible update and this updated 'z' variable is sent to the
            % members.
            if isempty(obj.coupled_variables)
                obj.coupled_variables = zeros(obj.n_members, size(val, 2));
            end
            obj.coupled_variables(obj.members_coupled(id), :) = val;
            obj.ctr = obj.ctr + 1;
            if obj.ctr == obj.n_members
                obj = obj.evalDualVariable();
            end
        end
        
        function [obj, message, epsilon] = evalCoupledVariables(obj, message)
            % Takes as input in 'message' the values of 
            % 
            % y + lambda/rho
            % 
            % of all of this Coordinators members. It returns in
            % 'message' the updates 'z' variable and in 'epsilon' a measure
            % the error in consensus of the 'z' variable.
            rows = size(message, 1);
            if rows == 0
                message = [];
                epsilon = 0;
                return;
            end
            
            if isempty(obj.z)
                z_prev = 0;
            else
                z_prev = obj.z;
            end
            
            for i=1:rows
                obj = obj.updateCoupledVariables(message(i, 1), message(i, 2:end));
            end
            
            message = [obj.members , (ones(size(obj.members))*obj.z)];
            
            temp = obj.z - z_prev;
            epsilon = norm(temp, 2);
        end
        
        function obj = update(obj)
            % Moves this Coordinator object to the next time-step.
            if ~isempty(obj.z)
                obj.z(1:(end-1)) = obj.z(2:end);
                obj.z(end) = 2*obj.z(end-1) - obj.z(end-2);
            end
            obj.k = obj.k + 1;
        end
        
        function obj = evalDualVariable(obj)
            % Updates the 'z' variable
            obj.z = sum(obj.coupled_variables, 1) ./ 2;
            obj.ctr = 0;
        end
    end
end
            