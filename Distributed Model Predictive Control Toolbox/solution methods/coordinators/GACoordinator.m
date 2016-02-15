classdef GACoordinator < Coordinator
    % A GACoordinator manages the communication of coupled varaible and dual
    % variable between two DualObjects. A GACoordinator object manages
    % the update of the dual variable, the update is performed according to
    % Nesterovs accelerated gradient method [1].
    % 
    % [1] - Y. Nesterov, "Introductory lectures on convex optimization",
    %       2004, Springer.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        % lambda and coupled variable properties
        ctr;
        lambda;
        alpha;
        mu;
        gamma;
        method; % Determines update rule
                % = 1 - fast gradient method
                % = 0 - gradient method
    end
    
    properties(Constant)
        FASTGRADIENT = 'fastgradient';
        GRADIENTASCENT = 'gradientascent';
    end
    
    methods
        function obj = GACoordinator(group) 
            % See help GACoordinator
            obj@Coordinator(group); %Call super class constructor
            obj.ctr = 0;
            obj.mu = 0;
            obj.gamma = 0;
            obj.k = 0;
            obj = obj.update();
            obj.method = 1;
        end
        
        
        function obj = setStepSize(obj, alpha)
            obj.alpha = alpha;
        end
        
        function obj = setMethod(obj, method)
            % Determines update rule
            % method = 1 - fast gradient method
            %        = 0 - gradient method
            if isnumeric(method) && method == 1
                obj.method = 1;
            elseif strcmp(method, GACoordinator.FASTGRADIENT)
                obj.method = 1;
            else
                obj.method = 0;
            end
        end
        
        function method = getMethod(obj)
            if ~isnumeric(obj.method)
                error('Method has not been set')
            end
            if obj.method == 1
                method = GACoordinator.FASTGRADIENT;
            else
                method = GACoordinator.GRADIENTASCENT;
            end
        end
        
        function obj = updateCoupledVariables(obj, id, val)
            % Stores the values of coupled variable for this Coordinators
            % members. When all the members has sent a new value of ther
            % coupled variable those values are used to compute the dual
            % varaible update and this updated dual variable is sent to the
            % members.
            if isempty(obj.coupled_variables)
                obj.coupled_variables = zeros(obj.n_members, size(val, 2));
                obj.lambda = zeros(1, size(val, 2));
            end
            obj.coupled_variables(obj.members_coupled(id), :) = val;
            obj.ctr = obj.ctr + 1;
            if obj.ctr == obj.n_members
                obj = obj.evalDualVariable();
            end
        end
        
        function [obj, message, epsilon] = evalCoupledVariables(obj, message)
            % Takes as input in 'message' the values of the coupled
            % variable of all of this Coordinators members. It returns in
            % 'message' the updates dual variable and in 'epsilon' a measure
            % the error in consensus of the coupled variable.
            rows = size(message, 1);
            if rows == 0
                message = [];
                epsilon = 0;
                return;
            end
            
            for i=1:rows
                obj = obj.updateCoupledVariables(message(i, 1), message(i, 2:end));
            end
            
            message = [obj.members , [obj.lambda; -obj.lambda]];
            
            epsilon = norm((obj.coupled_variables(1, :)-obj.coupled_variables(2, :)), 2);
            %disp(['Epsilon: ' num2str(epsilon)])
        end
        
        function obj = update(obj)
            % Moves this Coordinator object to the next time-step.
            if ~isempty(obj.lambda)
                obj.lambda(1:(end-1)) = obj.lambda(2:end);
                obj.lambda(end) = 2*obj.lambda(end-1) - obj.lambda(end-2);
            end
            obj.k=obj.k+1;
            obj.gamma = (sqrt(5)-1)/2;
        end
        
        function obj = evalDualVariable(obj)
            % Updates the dual variable
            if ~isnumeric(obj.alpha)
                error('The step size of the dual variable update has not been set')
            end
            if obj.method == 1
                obj = obj.fastGradient();
            else
                obj = obj.gradientAscent();
            end
            obj.ctr = 0;
        end
        
        function obj = fastGradient(obj)
            % Updates the dual variable with Nesterovs accelerated gradient
            % algorithm (see help FastGradient).
            mu1 = obj.lambda + ...
                (obj.coupled_variables(1, :)-obj.coupled_variables(2, :)) * obj.alpha;
            gamma1 = obj.gamma*(sqrt(obj.gamma*obj.gamma + 4) - obj.gamma) / 2;
            
            beta = obj.gamma*(1-obj.gamma) / (obj.gamma*obj.gamma + gamma1);
            
            obj.gamma = gamma1;
            
            obj.lambda = mu1 + beta*(mu1-obj.mu);
            obj.mu = mu1;
            
        end
       
        function obj = gradientAscent(obj)
            % Updates the dual variable using gradient ascent.
            obj.lambda = obj.lambda + ...
                (obj.coupled_variables(1, :)-obj.coupled_variables(2, :)) * obj.alpha;
        end
    end
end
            