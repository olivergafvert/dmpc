clear all; close all; clc;

% A test that uses the classes FastGradient, MultiStep and ADMM to 
% solve the distributed MPC problem of controlling five interconnected tanks.

%Setup system object with prediction horizon 20 and control horizon 10
sys = DmpcSys(20, 10);


A = 1;
B = [1 -1 -1 1];
C = 1;

LTI1 = ss(A, B, C, [], 1);

params.Q = 1;
params.R = 0.1*eye(4);
params.umax = [10; 10; 5; 2];
params.umin = [0; 0; -5; -2];
params.ymax = 10;
params.ymin = -10;
params.x_0 = -5;

params.soft = 3;

[sys, id1] = sys.addSubsystem(LTI1, params);

A = 1;
B = [1 -1 1];
C = 1;

LTI2 = ss(A, B, C, [], 1);

params.Q = 1;
params.R = 0.1*eye(3);
params.umax = [10; 2; 5];
params.umin = [0; -2; -5];
params.ymax = 10;
params.ymin = -1;
params.x_0 = 5;

[sys, id2] = sys.addSubsystem(LTI2, params);

A = 1;
B = [-1 1 1];
C = 1;

LTI3 = ss(A, B, C, [], 1);

params.Q = 1;
params.R = 0.1*eye(3);
params.umax = [2; 2; 5];
params.umin = [-2; -2; -5];
params.ymax = 10;
params.ymin = -10;
params.x_0 = 10;

[sys, id3] = sys.addSubsystem(LTI3, params);


A = 1;
B = [-1 1 -1];
C = 1;

LTI4 = ss(A, B, C, [], 1);

params.Q = 1;
params.R = 0.1*eye(3);
params.umax = [5; 5; 5];
params.umin = [-5; -5; -5];
params.ymax = 10;
params.ymin = -10;
params.x_0 = 10;

[sys, id4] = sys.addSubsystem(LTI4, params);


A = 1;
B = [1 -1];
C = 1;

LTI5 = ss(A, B, C, [], 1);

params.Q = 1;
params.R = 0.1*eye(2);
params.umax = [10; 5];
params.umin = [0; -5];
params.ymax = 100;
params.ymin = -100;
params.x_0 = -10;

[sys, id5] = sys.addSubsystem(LTI5, params);


% Connect the subsystems to each other

%sys = sys.connect(id1, 3, 'input', id2, 3, 'input');
sys = sys.connect(id1, 4, 'input', id3, 1, 'input');
sys = sys.connect(id2, 2, 'input', id3, 2, 'input');

%sys = sys.connect(id3, 3, 'input', id4, 1, 'input');
sys = sys.connect(id5, 2, 'input', id4, 2, 'input');
sys = sys.connect(id5, 1, 'input', id1, 2, 'input');
%sys = sys.connect(id2, 1, 'output', id3, 1, 'output');

sys.print()


% Try different solution methods

maxIter = 500;
thresh = 0.001;
n_simulations = 10;

iter = [];

disp('Fast Gradient Method')

sol = FastGradient(sys, maxIter, thresh, Solver.ADMM, Log.NONE);


sol = sol.sim(n_simulations);

iter = [iter ; sol.histIter];

disp('Multi-step gradient method')

sol = MultiStep(sys, maxIter, thresh, Solver.ADMM, Log.NONE);

sol = sol.sim(n_simulations);

iter = [iter ; sol.histIter];

disp('ADMM')

sol = ADMM(sys, maxIter, thresh, Solver.ADMM, Log.NONE);

sol = sol.sim(n_simulations);

iter = [iter ; sol.histIter];

plot(1:n_simulations, iter)

legend('Fast Gradient', 'Multi-step', 'ADMM')
