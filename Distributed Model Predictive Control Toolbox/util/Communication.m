classdef Communication
    % This object to simulates communication in a distributed mpc system.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        in; % pass message to a dmpc object
        out; % pass message to a group
    end
    
    methods
        function obj = Communication()
            % See help Communication
            obj.in = cell(1, 1);
            obj.out = cell(1, 1);
        end
        
        function obj = toGroup(obj, id, message)
            % id - the id of some mpc object
            % message - a matrix concisting of columns [id of dmpc ,
            % lambda]
            if size(message, 1) == 0
                return;
            end
            ids = message(:, 1);
            if isa(ids, 'cell')
                tmp = [];
                for i=1:length(ids)
                    tmp = [tmp; ids{i}];
                end
                ids = tmp;
            end
            for i=1:size(ids, 1)
                try
                    obj.in{ids(i)};
                catch why
                    obj.in{ids(i)} = [];
                end
                obj.in{ids(i)} = [obj.in{ids(i)} ; [id message(i, 2:end)]];
            end
        end
        
        function obj = toDMPC(obj, id, message)
            % id - the id of some mpc object
            % message - a matrix concisting of columns [id of group ,
            % lambda]
            ids = message(:, 1);
            if isa(ids, 'cell')
                tmp = [];
                for i=1:length(ids)
                    tmp = [tmp; ids{i}];
                end
                ids = tmp;
            end
            for i=1:size(ids, 1)
                try
                    obj.out{ids(i)};
                catch
                    obj.out{ids(i)} = [];
                end
                obj.out{ids(i)} = [obj.out{ids(i)} ; [id message(i, 2:end)]];
            end
        end
        
        function message = getMessageGroup(obj, id)
            try
                message = obj.in{id};
            catch
                message = [];
            end
        end
        
        function message = getMessageDMPC(obj, id)
            try
                message = obj.out{id};
            catch
                message = [];
            end
        end
        
        function obj = clearOut(obj)
            obj.out = cell(1, 1);
        end
        
        function obj = clearIn(obj)
            obj.in = cell(1, 1);
        end
    end
end
    