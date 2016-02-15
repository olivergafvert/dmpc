classdef DmpcSubsystem
    % A DmpcSubsystem object is a representation of an Mpc object which is
    % coupled to other Mpc objects within a distributed system. Each
    % DmpcSubsystem object contains an Mpc object and information about how
    % the Mpc object is coupled to other DmpcSubsystem objects. The
    % coupling of a variable between two DmpcSubsystem objects is 
    % represented by that the variable in each subsystem is members of the
    % same DmpcGroup object (see help DmpcGroup).
    
    %   Author: Oliver Gäfvert
    %   Email: oliverg@kth.se
    properties
        id; % The id of this DmpcSubsystem within the distributed system
        mpc; % The Mpc object
        coupled_variables; % a map connecting a coupled variable to its group
        groups; % A vector of group id's of the groups this subsystem belongs to
        n_groups; % The number of groups this subsystem belongs to
    end
    
    methods
        function obj = DmpcSubsystem(mpc)
            % See help DmpcSubsystem
            obj.id = ID.nextID();
            obj.mpc = mpc;
            obj.groups = zeros(1, 1);
            obj.n_groups = 0;
            obj.coupled_variables = zeros(size(mpc.LTI.c, 1), 1);
        end
        
        function obj = addGroup(obj, groupID, coupled_variableID)
            % Adds the group with group id 'groupID' to this DmpcSubsystem
            % object. The group id is added to the property 'groups' and
            % the 'coupled_variableID' is added to the property
            % coupled_variables. If the property coupled_variables already
            % contatins a value at position 'coupled_varialbeID' a new
            % column of zeros is created and the value 'groupID' is added
            % at position 'coupled_varialbeID' in this new column. Hence if
            % the system uses node varialbes each row in the matrix
            % property coupled_variables contains the groups each output
            % variable is coupled to. If the system uses edge varialbes the
            % property coupled_variables will be a vector.
            obj.n_groups = obj.n_groups +1;
            obj.groups(obj.n_groups) = groupID;
            
            found = 0;
            if coupled_variableID > size(obj.coupled_variables, 1)
                obj.coupled_variables = [obj.coupled_variables; zeros(...
                    coupled_variableID-size(obj.coupled_variables, 1), size(obj.coupled_variables, 2))];
            end
            for i=1:size(obj.coupled_variables, 2)
                if obj.coupled_variables(coupled_variableID, i) == 0
                    obj.coupled_variables(coupled_variableID, i) = groupID;
                    found = 1;
                    break;
                end
            end
            
            if found == 0
                obj.coupled_variables = [obj.coupled_variables, zeros(size(obj.coupled_variables, 1), 1)];
                obj.coupled_variables(coupled_variableID, end) = groupID;
            end
        end
        
        function row = getCoupledVar(obj, groupID)
            % Returns the index of the output variable that is coupled to
            % the group with group id 'groupID'.
            [row, col] = find(obj.coupled_variables == groupID);
        end
        
        function coupled = isCoupled(obj, y_index)
            % Returns 1 if the output variable with index 'y_index' is
            % coupled to some group, it returns 0 else.
            if obj.coupled_variables(y_index, 1) == 0
                coupled = 1;
            else
                coupled = 0;
            end
        end
        
        function print(obj)
            % Prints some information about this DmpcSubsystem object
            disp(['Subsystem id: ' num2str(obj.id)]);
            disp('Coupled variables: ');
            disp(obj.coupled_variables);
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
        
        function obj = clearGroups(obj)
            % Clears the couplings of the DmpcSubsystem object.
            obj.n_groups = 0;
            obj.groups = cell(1, 1);
            obj.coupled_variables = [];
        end
    end
end
            