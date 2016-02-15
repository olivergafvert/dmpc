classdef Mpc
    % Initiates an Mpc object with fields
    % mpc.LTI - an LTI object
    %    .N - prediction horizon
    %    .Nc - control horizon
    %    .Q - weight matrix for x
    %    .R - weight matrix for u
    %    .x_k - current x-value
    %    .umax - maximum allowed u-value(s)
    %    .umin - minimum allowed u-value(s)
    %    .ymax - maximum allowed y-value(s)
    %    .ymin - minimum allowed y-value(s)
    %    .H - Hessian in the qp-problem solved by Mpc
    %    .F - the matrix in linear term in qp-problem solved by Mpc
    %    .Ac - constraint matrix in qp-problem solved by Mpc
    %    .b - constraint matrix in qp-problem solved by Mpc
    %    .Bx - constraint matrix in qp-problem solved by Mpc
    %    .Conv - Convolusion matrix (see help Mpc.initWeights)
    %    .M - matrix with exponents of LTI.a (see help Mpc.initWeights)
    %    .soft -  determines soft constrsints
    %             = 0 - no soft constraints
    %             = 1 - soft output constrints
    %             = 2 - soft input constraints
    %             = 3 - both soft input and output constraints
    %    
    % and computes the matrices needed to solve the quadratic programming
    % problem
    % 
    % min (1/2)*u'*H*u + (F*mpc_k)'*u
    %  x
    % s.t Ac*u <= b + Bx*x_k
    % 
    % (see help Mpc.initWeights for more detailed information)
    % 
    % To create an Mpc object the following need to be specified:
    % 
    % Mpc(LTI, params)
    % 
    % LTI - an LTI object
    % params - an object with parameters. It should contain the fields
    %          N, Nc, Q, R, x_0, umax, umin, ymax, ymin, soft (optinal).
    % 
    % Given an LTI object and an 'params' struct an Mpc object is created
    % as follows
    % 
    % obj = Mpc(LTI, params)
    % 
    % the input signals for the next time step can then be computed by e.g
    % Matlabs quadprog solver
    % 
    % u = quadprog(obj.H, obj.F*obj.x_k, obj.Ac, obj.b+obj.Bx*obj.x_k)
    % 
    % note that if soft constraints are used the 'u'-vector also contains
    % slack variables. to remove these slack variables set
    % 
    % bc = size(obj.LTI.b, 2)
    % u = u(1:bc*obj.Nc)
    
    %   Author: Oliver Gäfvert
    %   E-mail: oliverg@kth.se
    properties
        LTI; % The LTI object of the system
        N; % prediection horizon
        Nc; % control horizon
        Q; % Weight matrix in cost function
        R; % Weight matrix in cost function
        x_k; % current state of the system
        umax; % maximum constraints on u (input signal)
        umin; % minimum constraints on u (input signal)
        ymax; % maximum constraints on y (output signal)
        ymin; % minimum constraints on y (output signal)
        H; % Hessian in the qp-problem
        F; % the matrix in linear term in qp-problem
        Ac; % constraint matrix
        b; % constraint matrix
        Bx; % constraint matrix
        Conv; % Convolusion matrix
        M; % M-matrix (the matrix with expoents of A=LTI.a) (C=LTI.c)
        soft; % information about if soft constraints are used
    end
    
    methods
        function mpc = Mpc(LTI, params)
            % See help Mpc
            mpc.LTI = LTI;
            mpc.N = params.N;
            mpc.Nc = params.Nc;
            mpc.Q = params.Q;
            mpc.R = params.R;
            mpc.x_k = params.x_0;
            mpc.umax = params.umax;
            mpc.umin = params.umin;
            mpc.ymax = params.ymax;
            mpc.ymin = params.ymin;
            if isfield(params, 'soft')
                mpc = mpc.initWeights(params.soft);
            else
                mpc = mpc.initWeights(0);
            end
        end
        
        function mpc = initWeights(mpc, soft)
            % Incorperates the weight matrices and constraints needed for
            % optimization into the Mpc object.
            % 
            % Input
            % soft - determines soft constrsints
            %               = 0 - no soft constraints
            %               = 1 - soft output constrints
            %               = 2 - soft input constraints
            %               = 3 - both soft input and output constraints
            % 
            % Weight matrices for objective function (1/2)*u'*H*u + (F*mpc_k)'*u
            % mpc.H - H
            % mpc.F - F 
            % 
            % H is computed as H = C'*\tilde(Q)*C + \tilde(R) where C is the
            % convolusional matrix, \tilde(Q) and \tilde(R) are large block
            % matrices with mpc.Q and mpc.R in the diagonal respectively.
            % F is computed as F = C'*\tilde(Q)*M where M is a matrix concisting of
            % exponants of mpc.LTI.a.
            % 
            % When soft constraints are used slack variables are added to the
            % optimization problem. H and F are modified (into Hs and Fs) by appending
            % matricies to support slack variables. The slack variables are integrated
            % into u so that u = [u ; e] where e is a vector of slack variables.
            % 
            % Linear constraints Ac*u < b + Bx*x_k. The constraints mpc.umax, mpc.umin,
            % mpc.xmax, mpc.xmin are reformulated into this matrix form.
            % mpc.Ac - Ac
            % mpc.b - b
            % mpc.Bx - Bx
            % 
            % If soft constrints are used the additional constraint e >= 0 is imposed 
            % on the slack variables e. These constraints are included in the matrices
            % Ac, b and Bx.

            % Init convolusion matrix, C
            [br, bc] = size(mpc.LTI.b);
            C = zeros(mpc.N*br, mpc.Nc*bc);
            C(1:br, 1:bc) = mpc.LTI.b;

            % Init matrix M concisting of exponants of mpc.LTI.a
            [ar, ac] = size(mpc.LTI.a);
            M = [mpc.LTI.a ; zeros((mpc.N-1)*ar, ac)]; 
            tmp_A = mpc.LTI.a;

            R = mpc.R;

            % Computing H and F
            H = 0;
            F = 0;

            %if mpc.N > mpc.Nc
                K = -dlqr(mpc.LTI.a, mpc.LTI.b, mpc.Q, mpc.R); %LQ-optimal gain matrix
            

            A_BK = mpc.LTI.a+mpc.LTI.b*K;
            tmp_A_BK = A_BK;

            for i=1:(mpc.N-1)
                rows = (i*ar+1):((i+1)*ar);
                if i < mpc.Nc
                    R = blkdiag(R, mpc.R); %create large matrix with mpc.R in the diagonal
                    C(rows, :) = [tmp_A*mpc.LTI.b, C(rows-br, 1:end-bc)];
                    tmp_A = mpc.LTI.a*tmp_A;
                    M(rows, :) = tmp_A;
                else
                    C(rows, :) = tmp_A_BK*C(((mpc.Nc-1)*ar+1):(mpc.Nc*ar), :);
                    M(rows, :) = tmp_A_BK*tmp_A;
                    tmp_A_BK = A_BK*tmp_A_BK;
                end

                rows = rows - ar;

                CQ = C(rows, :)'*mpc.Q;
                H = H + CQ*C(rows, :);
                F = F + CQ*M(rows, :); 
            end


            % Computing the terminal weight matrix for use in H and F

            Qt = dlyap(A_BK',mpc.Q+K'*mpc.R*K); %calculate terminal weight matrix

            rows = ((mpc.N-1)*ar+1):(mpc.N*ar);

            CQ = C(rows, :)'*Qt;
            H = H + CQ*C(rows, :) + R;
            F = F + CQ*M(rows, :); 


            mpc.H = H;
            mpc.F = F;

            % mpc.H
            eps = 2.22e-16;
            if norm(mpc.H-mpc.H',inf) > eps
                disp(['H is not symmetric, norm(mpc.H-mpc.H^T,inf) = ', num2str(norm(mpc.H-mpc.H',inf))])
                disp('Resetting H = (H^T + H)/2');
                mpc.H = (mpc.H+mpc.H')/2;
            end


            % Computing Ac and b

            [cr, cc] = size(mpc.LTI.c);

            umax = ones(mpc.Nc*bc, 1); umin = ones(mpc.Nc*bc, 1);
            ymax = ones(mpc.N*cr, 1); ymin = ones(mpc.N*cr, 1);

            % K*C*u <= umax - K*M*x_k
            uTmax = zeros((mpc.N-mpc.Nc)*bc, 1);
            uTmin = zeros((mpc.N-mpc.Nc)*bc, 1);
            KC = zeros(bc*(mpc.N-mpc.Nc), size(C, 2));
            KM = zeros(bc*(mpc.N-mpc.Nc), ac);

            CC = [];

            for i=1:mpc.N
                CC = blkdiag(CC, mpc.LTI.c);
                if i <= mpc.Nc
                    rows = ((i-1)*bc+1):(i*bc);
                    umax(rows) = mpc.umax;
                    umin(rows) = mpc.umin;
                else
                    crows = ((i-1)*ar+1):(i*ar);
                    rows = ((i-mpc.Nc-1)*bc+1):((i-mpc.Nc)*bc);
                    uTmax(rows) = mpc.umax;
                    uTmin(rows) = mpc.umin;
                    KM(rows, :) = K*M(crows, :);
                    KC(rows, :) = K*C(crows, :);
                end

                rows = ((i-1)*cr+1):(i*cr);
                ymax(rows) = mpc.ymax;
                ymin(rows) = mpc.ymin;
            end

            % Compute the D-matrix
            DD = mpc.LTI.d;
            for i=1:(mpc.Nc-1)
                DD = blkdiag(DD, mpc.LTI.d);
            end
            
            
            mpc.b = [umax ; -umin ; ymax ; -ymin ; uTmax ; -uTmin];

            Y = CC*M;

            mpc.Bx = [zeros(mpc.Nc*bc*2, ac) ; -Y ; Y ; -KM ; KM];

            
            % (C*Conv + DD)*u == y
            CC_C = CC*C;

            mpc.Conv = C(1:(br*mpc.Nc), :);
            mpc.M = M(1:(br*mpc.Nc), :);

            mpc.Ac = [eye(mpc.Nc*bc) ; -eye(mpc.Nc*bc) ; CC_C ; -CC_C ; KC ; -KC];

          
            if soft >= 1
                mpc.soft = soft;
                if soft >= 2 % init slack variables in input signals
                    E = eye(bc); %using identity matrix as weight matrix for slack variables
                    EE = [];
                    for i=1:mpc.N
                        EE = [EE ; E];
                    end
                    Ec = [-EE(1:(mpc.Nc*bc), :) ; -EE(1:(mpc.Nc*bc), :) ; zeros(2*size(CC_C, 1), size(E, 2)) ; -EE(1:((mpc.N-mpc.Nc)*bc), :) ; -EE(1:((mpc.N-mpc.Nc)*bc), :)];
                end

                if soft ~= 2 % init slack vairables in output signals
                    V = eye(cr); %using identity matrix as weight matrix for slack variables
                    VV = [];
                    for i=1:mpc.N
                        VV = [VV ; V];
                    end
                    Vc = [zeros(mpc.Nc*bc*2, size(VV, 2)) ; -VV ; -VV ; zeros(size(KC, 1)*2, size(VV, 2))];
                end 

                % Modify H, F, Ac, b and Bx to support slack variables
                if soft > 2
                    Acs = [[mpc.Ac , Vc, Ec] ; [zeros(cr, size(mpc.Ac, 2)) , -eye(cr), zeros(cr, size(E, 2))] ; [zeros(bc, size(mpc.Ac, 2)+size(V, 2)) , -eye(bc)]]; 
                    bs = [mpc.b ; zeros(cr, 1) ; zeros(bc, 1)];
                    Bxs = [mpc.Bx ; zeros(cr, size(mpc.Bx, 2)) ; zeros(bc, size(mpc.Bx, 2))];
                    Hs = blkdiag(blkdiag(H, V), E);
                    Fs = [mpc.F ; zeros(size(V, 1), size(mpc.F, 2)) ; zeros(size(E, 1), size(mpc.F, 2))];
                elseif params.soft == 2
                    Acs = [[mpc.Ac , Ec] ; [zeros(bc, size(mpc.Ac, 2)) , -eye(bc)]]; 
                    bs = [mpc.b ; zeros(bc, 1)];
                    Bxs = [mpc.Bx ; zeros(bc, size(mpc.Bx, 2))];
                    Hs = blkdiag(H, E);
                    Fs = [mpc.F ; zeros(size(E, 1), size(mpc.F, 2))];
                else
                    Acs = [[mpc.Ac , Vc] ; [zeros(cr, size(mpc.Ac, 2)) , -eye(cr)]]; 
                    bs = [mpc.b ; zeros(cr, 1)];
                    Bxs = [mpc.Bx ; zeros(cr, size(mpc.Bx, 2))];
                    Hs = blkdiag(H, V);
                    Fs = [mpc.F ; zeros(size(V, 1), size(mpc.F, 2))];
                end
                mpc.H = Hs;
                mpc.F = Fs;
                mpc.Ac = Acs;
                mpc.b = bs;
                mpc.Bx = Bxs;
            else
                mpc.soft = 0;
            end
        end
        
        function obj = reCompile(obj)
            % Recompiles the weights of the Mpc object (should be done if 
            % paramters are changed etc...).
            obj = obj.initWeights(obj.soft);
        end
    end
end


        
        
        