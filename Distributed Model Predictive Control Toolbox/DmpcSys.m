classdef DmpcSys
    % The DmpcSys class represents a distributed system of Mpc controllers.
    % It represents the subsystems as DmpcSubsystem objects and coupling
    % between subsystems by DmpcGroup objects.
    % 
    % To create a DmpcSys object the following need to be specified:
    % 
    % 'N' - Prediction horizon
    % 'Nc' - control horizon 
    % 'type' - the type of architecture of the system (if it
    %          is modeled by edge variables or node variables). 
    %          'type' can take the values 'Edge' for edge variables
    %          and 'Node' for node varialbes. When the system is
    %          modeled by edge variables each output varialbe coupling
    %          can only be coupled between two subsystems. If the coupling
    %          is between more than two subsystems a 'virtual'
    %          output is created by copying the coupled output 
    %          variable into a new output variable. When using node
    %          varialbes one output variable can be coupled to
    %          several subsystems.
    % 
    % The following command will create a DmpcSys object that has N = 10,
    % Nc = 5 and uses node variables in the coupled variables.
    % 
    % sys = DmpcSys(10, 5, 'Node')
    
    %   Author: Oliver Gäfvert
    %   Email: oliverg@kth.se
    properties
        subsID_subs; % a map between subsystem ids and subsystem positions in the subsystems array
        groupID_group; % a map between group ids and positions of groups in the groups cell array
        subsystems; % a cell array containing the DmpcSubsystem objects
        n_subsystems; % the number of subsystems in the system
        groups; % a cell array containing the DmpcGroup objects
        n_groups; % the number of groups in the system
        N; %prediction horizon of the system
        Nc; %control horizon of the system
        c_type; %Connection type (1 - Edge variables, 0 - Node variables)
    end
    
    methods
        function obj = DmpcSys(N, Nc, type)
            % See help DmpcSys
            obj.N = N;
            obj.Nc = Nc;
            obj.subsID_subs = zeros(1, 1);
            obj.groupID_group = zeros(1, 1);
            obj.subsystems = cell(1, 1);
            obj.n_subsystems = 0;
            obj.groups = cell(1, 1);
            obj.n_groups = 0;
            if nargin > 2
                if strcmp(type, 'Node')
                    obj.c_type = 0;
                else
                    obj.c_type = 1;
                end
            else
                obj.c_type = 1;
            end   
        end
        
        function pos = getGroupPosition(obj, id)
            % Returns the position of the group with group id 'id' in the
            % cell array property 'groups'.
            pos = obj.groupID_group(id);
        end
        
        function pos = getSubsystemPosition(obj, id)
            % Returns the position of the subsystem with id 'id' in the
            % cell array property 'subsystems'.
            pos = obj.subsID_subs(id);
        end
        
        function [obj, id] = addSubsystem(obj, mpc, params)
            % [obj, id] = obj.addSubsystem(LTI, params) adds the subsystem
            % to this DmpcSys object by fist creating an Mpc object using
            % the parameters 'LTI' and 'params'. See help Mpc for what
            % should be contained in the 'params' struct. When a Mpc object
            % has been created a DmpcSubsystem object is created given
            % the Mpc object as a parameter and added to the DmpcSys object.
            % 
            % [obj, id] = obj.addSubsystem(mpc) creates a DmpcSubsystem 
            % object given the parameter mpc (which should be an Mpc
            % object) and adds this to the DmpcSys object.
            if nargin > 2
                params.N = obj.N;
                params.Nc = obj.Nc;
                mpc = Mpc(mpc, params);
            end
            if mpc.Nc ~= obj.Nc
                mpc.Nc = obj.Nc;
                mpc.N = obj.N;
                mpc = mpc.reCompile();
            end
            subsys = DmpcSubsystem(mpc);
            id = subsys.id;
            obj.n_subsystems = obj.n_subsystems + 1;
            obj.subsystems{obj.n_subsystems} = subsys;
            obj.subsID_subs(subsys.id) = obj.n_subsystems; % add subsystem id to map
        end
       
        function obj = addGroup(obj, subsIDs, subsIDYs)
            % A group is a virtual link between two (or more) subsystems representing a
            % coupled variable.
            group = DmpcGroup(subsIDs, subsIDYs); % create a group for the coupled variable
            
            % add this group to all subsystems belonging to the group
            for i=1:max(size(subsIDs))
                pos = obj.subsID_subs(subsIDs(i));
                obj.subsystems{pos} = obj.subsystems{pos}.addGroup(group.id, subsIDYs(i));
            end
            
            obj.n_groups = obj.n_groups + 1;
            obj.groups{obj.n_groups} = group;
            obj.groupID_group(group.id) = obj.n_groups;
        end
        
        function obj = connect(obj, s1ID, s1y, type1, s2ID, s2y, type2)
            % Connects two subsystems via output or input signals. Coupling
            % in input signals is represented internally as coupling in
            % output signals This is done by adding output signals to the 
            % LTI models of the subsytems encoding this coupling.
            %
            % The coupling of two output (or input) variables is
            % represented internally as the variables being a member of the
            % same DmpcGroup object (see help DmpcGroup).
            % 
            % s1ID, S2ID - could be either a DmpcSubsystem object or the id of a
            %              DmpcSubsystem object.
            % s1y, s2y - the index of the output or input signal in
            %            subsystem (with id) s1ID, s2ID.
            % type1, type2 - can take the values 'input' for a coupling in
            %                the input signals and 'output' for coupling in
            %                output signals. (Default is 'output')
            %
            % obj = obj.connect(s1ID, s1y, s2ID, s2y) connects the output
            % signal with index 's1y' in the subsystem with id 's1ID' to
            % the output signal with index 's2y' in the subsystem with id
            % 's2ID'.
            % 
            % obj = obj.connect(s1ID, s1y, type1, s2ID, s2y, type2) connects 
            % the output/input signal (depending on if type1='input' or 
            % type1='output') with index 's1y' in the subsystem with id 
            % 's1ID' to the output/input signal (depending on if 
            % type1='input' or type1='output') with index 's2y' in the 
            % subsystem with id 's2ID'.
            if nargin < 6
                s2y = s2ID;
                s2ID = type1;
                type1 = 'output';
                type2 = type1;
            end
            
            if isa(s1ID, 'DmpcSubsystem')
                s1ID = s1ID.id;
            end
            if isa(s2ID, 'DmpcSubsystem')
                s2ID = s2ID.id;
            end
            
            y_index = [s1y; s2y];
            type = cell(2, 1);
            type{1} = type1;
            type{2} = type2;
            id = [s1ID; s2ID];
            yid = [s1y s2y];
            
            for i = 1:max(size(id))
                if strcmp(type{i}, 'input')
                    pos = obj.getSubsystemPosition(id(i));
                    
                    Ts = obj.subsystems{pos}.mpc.LTI.Ts;
                    
                    A = obj.subsystems{pos}.mpc.LTI.a;
                    B = obj.subsystems{pos}.mpc.LTI.b;
                    
                    % Resize C to create a new y
                    C = [obj.subsystems{pos}.mpc.LTI.c ; ...
                        zeros(1, size(obj.subsystems{pos}.mpc.LTI.c, 2))];
                    y_index(i) = size(C, 1);

                    % Resize D to make y_new = u(s1y)
                    D = [obj.subsystems{pos}.mpc.LTI.d ; ...
                         zeros(1, size(obj.subsystems{pos}.mpc.LTI.d, 2))];
                  
                    obj.subsystems{pos}.mpc.LTI = ss(A, B, C, D, Ts); 
                    obj.subsystems{pos}.mpc.LTI.d(end, yid(i)) = 1;
                
                elseif strcmp(type{i}, 'output')
                    pos = obj.getSubsystemPosition(id(i));
                    if obj.subsystems{pos}.isCoupled(y_index(i)) && obj.c_type == 1 % Only add new output if c_type = 1 (edge variables)
                        Ts = obj.subsystems{pos}.mpc.LTI.Ts;
                    
                        A = obj.subsystems{pos}.mpc.LTI.a;
                        B = obj.subsystems{pos}.mpc.LTI.b;

                        % Resize C to create a new y
                        C = obj.subsystems{pos}.mpc.LTI.c;
                        C = [C;C(y_index(i), :)];
                        
                        % Resize D to make y_new = u(s1y)
                        D = obj.subsystems{pos}.mpc.LTI.d;
                        D = [D;D(y_index(i), :)];
                        
                        y_index(i) = size(C, 1);

                        obj.subsystems{pos}.mpc.LTI = ss(A, B, C, D, Ts); 
                    end
                        
                else
                    disp(['Unknown type: ' type{i}]);
                end
            end
            
            obj = obj.addGroup(id, y_index);
            
        end
        
        function print(obj)
            % print(obj)
            % Prints information about the subsystems and connections contained in
            % the system.
            for i=1:obj.n_subsystems
                obj.subsystems{i}.print();
            end
            for i=1:obj.n_groups
                obj.groups{i}.print();
            end
        end
        
       
    end
end