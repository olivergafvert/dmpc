classdef Log
    % This class is used to log information in the Distributed Model Predictive
    % Control Toolbox.

    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        filename;
        buffer;
        n_buffer; %the number of elements in the buffer
        buffer_max_size;
        type;
    end
    
    properties(Constant)
        NONE = 'none';
        ALL = 'all';
        SUBSYS = 'subsys';
        COORDINATORS = 'coordinators';
        SOLVER = 'solver';
        EXTENSION = '.log';
        DIRECTORY = './';
    end
    
    methods
        function obj = Log(filename, type, buffer_size)
            obj.filename = filename;
            obj.n_buffer = 0;
            if nargin < 3
                obj.buffer_max_size = 100;
            else
                obj.buffer_max_size = buffer_size;
            end
            obj.buffer = cell(1, obj.buffer_max_size);
            obj.type = type;
            fileID = fopen([Log.DIRECTORY obj.filename Log.EXTENSION],'w');
            fprintf(fileID,'');
            fclose(fileID);
        end
        
        function obj = log(obj, data)
            obj.n_buffer = obj.n_buffer +1;
            obj.buffer{obj.n_buffer} = data;
            if obj.n_buffer == obj.buffer_max_size
                obj = obj.saveToFile();
            end
        end
        
        function obj = saveToFile(obj)
            fileID = fopen([Log.DIRECTORY obj.filename Log.EXTENSION],'a');
            for i=1:obj.n_buffer
                if ~isempty(obj.buffer{i})
                    fprintf(fileID,obj.buffer{i});
                end
            end
            fclose(fileID);
            obj.n_buffer = 0;
            obj.buffer = cell(1, obj.buffer_max_size);
        end
    end
end