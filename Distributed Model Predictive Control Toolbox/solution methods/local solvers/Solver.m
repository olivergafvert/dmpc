classdef Solver
    % A basic quadratic programming solver that uses Matlabs quadprog
    % solver.
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        H;
        f;
        A;
        b;
        Bx;
        c;
    end
    
    properties(Constant)
        QUADPROG = 'quadprog';
        ADMM = 'admm';
        ADMM_LOGGED = 'admm_logged';
    end
    
    methods
        function obj = Solver(H, f, A, b, Bx, x_k)
            obj.H = H;
            obj.f = f;
            obj.A = A;
            obj.b = b;
            obj.Bx = Bx;
            
            obj.c = obj.b + obj.Bx*x_k;
        end
        
        function obj = update(obj, x_k)
            obj.c = obj.b + obj.Bx*x_k;
        end
        
        function [obj, x, flag] = solve(obj, f, x)
            options = optimset('Display', 'off');

            n = nargin;
            if n < 3
                x = [];
                if n < 2
                    f = obj.f;
                end
                
                
            end
            
            warning('off', 'all')
            [x, fval, flag] = quadprog(obj.H, f, obj.A, obj.c, [], [], [], [], x, options);
            warning('on', 'all')
        end
    end
end