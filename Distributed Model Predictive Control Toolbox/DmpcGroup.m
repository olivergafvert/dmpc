classdef DmpcGroup
    % A DmpcGroup object is a representation of a coupled variable between
    % two DmpcSubsystem objects. If a variable in a DmpcSubsystem object is
    % coupled to another DmpcSubsystem object this coupling is represented
    % by the variable (in each subsystem) being a member of the same DmpcGroup
    % object.
    
    %   Author: Oliver Gäfvert
    %   Email: oliverg@kth.se
    properties
        id; % The id of this group within the distributed system
        members; % array of dmpc object ids
        n_members; % the number of members of this group
        coupledVarIds; % an array containing an identifier of the coupled variable in each member
    end
    
    methods
        function obj = DmpcGroup(members, coupledVarIds)
            % See help DmpcGroup
            obj.id = ID.nextID();
            obj.members = members;
            obj.n_members = max(size(members));
            obj.coupledVarIds = coupledVarIds;
        end
        
        function members = getMembers(obj)
            % Returns the members of this group.
            members = obj.members;
        end
        
        function print(obj)
            % Prints some information about this group.
            disp(['Group id: ' num2str(obj.id)])
            disp('Members: ')
            disp(obj.members)
            disp('Coupled variable identifier for each member: ')
            disp(obj.coupledVarIds)
        end
       
    end
end
            