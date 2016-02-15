classdef ID
    % This class is used to generate IDs. To following line generates a new
    % id
    % 
    % id = ID.nextID();
    % 

    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    methods( Static )
        function ret = nextID()
            persistent id;
            if isempty( id )
                id = 0;
            end
            id = id + 1;
            ret = id;
        end
    end
end
            