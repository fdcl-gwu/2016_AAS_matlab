% 19 July 2016
% run asteroid shooting on HPC
fprintf('STARTING ASTEROID SHOOTING\n')
start_time = tic;

%% ASTEROID CONSTANTS
constants = load_constants('castalia','true'); % only 1024 faces
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

constants.pot_model = 'polyhedron'; % or mascon or matlab

constants.ode_options = odeset('RelTol',1e-9,'AbsTol',1e-9);
constants.optfsolve = optimoptions(@fsolve,'Display','off','TolFun',1e-6,'TolX',1e-9,...
    'MaxIter',5000,'MaxFunEvals',5000, 'Algorithm', 'trust-region-reflective','Jacobian','off',...
    'DerivativeCheck','off', 'FunValCheck','on');

num_steps = 1000;
constants.num_steps = num_steps;
constants.num_seg = 2; % number of segments for multiple shooting

% DEFINE CENTER OF ROTATION FOR POINCARE CUT
% angles to define the poincare section
constants.center_vec = zeros(3,1);
constants.pmap = [1 3 4 6]; % Indices of Poincare map

% THINGS TO CHANGE
constants.um = 1e-7; % 0.4e-7 to 12e-7 for 20-600 mN electric thrusters on 500 kg spacecraft
constants.num_angles = 10; % discretization of Poincare section

constants.filename = 'hpc_stage1_hightol.mat';
%% INITIALIZE TRANSFER
% this generates the initial and desired orbits (not any of the
% subintervals)
fprintf('INITIALIZE TRANSFER\n');

% constants.control_switch = 'off';
% [t_initial,state_initial,t_target,state_target]=initialize_transfer(constants);

% initial step beginning
load('./results/initial_transfer4.mat');

tf = t_initial(end);
xcf = state_initial(end,1:6)';

% INITIAL GUESS OF COSTATE AND LAGRANGE MULTIPLIERS
x0_i = state_initial(end,:)';
h0_i = 1e-3*ones(6,1);
beta_i = 1e-3*ones(4,1);


%% CALL THE SHOOTING METHOD

[sol_output] = asteroid_shooting(x0_i,h0_i,beta_i,xcf,tf,constants);

end_time = toc(start_time);
fprintf('FINISHED ASTEROID SHOOTING IN %5.4f sec\n',end_time);
