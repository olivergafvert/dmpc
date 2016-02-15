classdef Coordinator
    % A Coordinator manages the communication of coupled varaibles and/or dual
    % variables between two DualObjects. A Coordinator object simply passes the
    % message sent by one DualObject to the DualObject coupled to that
    % object by this Coordinator object.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        id; % the id of this coordinator within the distributed system
        members; % vector of DualObject (or subclass) object ids
        n_members; % the number of subsystems that are coupled with this coordinator
        
        coupled_variables; % a matrix consisting of the value of the coupled
                           % varialbe in each subsystem (that is coupled by
                           % this variable).
        members_coupled; % a map mapping member id's to positions in coupled_varialbes.
        logger;
        k; % current time-step
    end
    
    methods
        function obj = Coordinator(group)
            % See help Coordinator
            obj.id = group.id;
            obj.members = group.members;
            obj.n_members = group.n_members;
            obj.members_coupled = zeros(1, 1);
            for i=1:obj.n_members
                obj.members_coupled(obj.members(i)) = i;
            end
            if obj.n_members ~= 2
                error('A Coordinator object should have exactly 2 members')
            end
            obj.logger = 0;
        end
        
        function obj = updateCoupledVariables(obj, id, val)
            % Stores the values of coupled variable for this Coordinators
            % members.
            if isempty(obj.coupled_variables)
                obj.coupled_variables = zeros(obj.n_members, size(val, 2));
            end
            obj.coupled_variables(obj.members(id), :) = val;
        end
        
        function [obj, message, epsilon] = evalCoupledVariables(obj, message)
            % Takes as input in 'message' the values of the coupled
            % variable of this Coordinators members. It swaps the
            % values of the coupled variables returns this in 'message',
            % the variable 'epsilon' is not used and therefore set to zero.
            epsilon = 0;
            rows = size(message, 1);
            if rows == 0
                message = [];
                return;
            end
            
            if isa(message(1, 1), 'cell')
                m1 = message(1, 1);
                m1 = m1{1};
                m2 = message(2, 1);
                m2 = m2{1};
            else
                m1 = message(1, 1);
                m2 = message(2, 1);
            end
            message = [[m1, message(2, 2:end)]; ...
                [m2, message(1, 2:end)]];
        end
        
        function obj = initLog(obj, type, ext)
            if nargin < 3
                ext = '';
            end
            if strcmp(type, Log.ALL) || strcmp(type, Log.COORDINATORS)
                obj.logger = Log([ext 'Coordinator' num2str(obj.id)], type);
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
        end
        
        function obj = update(obj)
            obj.k = obj.k + 1;
        end
    end
end
            