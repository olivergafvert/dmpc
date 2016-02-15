classdef MSCoordinator < GACoordinator
    % A Coordinator that handles the communication and dual variable update
    % between two subsystems using the Multi-Step Gradient Method.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        lambda_prev;
        beta;
    end
    
    methods
        function obj = MSCoordinator(group) 
            obj@GACoordinator(group); %Call super class constructor
            obj.lambda_prev = 0;
        end
        
        function obj = evalDualVariable(obj)
            % Updates the dual variable using a multi-step gradient method
            % so that 
            % 
            % lambda = lambda + alpha*(d(lambda)) + beta*(lambda -
            % lambda_prev)
            % 
            % where d(lambda) denotes the dual gradient and alpha and beta
            % are tunale parameters.
            prev = obj.lambda;
            obj.lambda = obj.lambda + ...
                (obj.coupled_variables(1, :) - obj.coupled_variables(2, :)) * obj.alpha ...
                + obj.beta*(obj.lambda-obj.lambda_prev);
            obj.lambda_prev = prev;
            obj.ctr = 0;
        end
        
        function obj = setStepSize(obj, alpha, beta)
            % Sets the alpha and beta that are used in the update rule (see
            % help MSCoordinator.evalDualVariable).
            obj.alpha = alpha;
            obj.beta = beta;
        end
        
        function obj = update(obj)
            % Moves this Coordinator object to the next time-step.
            if ~isempty(obj.lambda)
                obj.lambda(1:(end-1)) = obj.lambda(2:end);
                obj.lambda(end) = 2*obj.lambda(end-1) - obj.lambda(end-2);
            end
            if ~isempty(obj.lambda_prev)
                obj.lambda_prev(1:(end-1)) = obj.lambda_prev(2:end);
                obj.lambda_prev(end) = 2*obj.lambda_prev(end-1) - obj.lambda_prev(end-2);
            end
            obj.k=obj.k+1;
        end
    end
end
            