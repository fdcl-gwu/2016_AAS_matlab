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
constants.optfsolve = optimoptions(@fsolve,'Display','iter','TolFun',1e-6,'TolX',1e-9,...
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

constants.filename = 'hpc_stage3_z_hightol.mat';
%% INITIALIZE TRANSFER
% this generates the initial and desired orbits (not any of the
% subintervals)
fprintf('INITIALIZE TRANSFER\n');

% constants.control_switch = 'off';
% [t_initial,state_initial,t_target,state_target]=initialize_transfer(constants);

% % initial step beginning
% load('./results/initial_transfer5.mat');
% 
% tf = t_initial(end);
% xcf = state_initial(end,1:6)';
% 
% % INITIAL GUESS OF COSTATE AND LAGRANGE MULTIPLIERS
% x0_i = state_initial(end,:)';
% h0_i = 1e-3*ones(6,1);
% beta_i = 1e-3*ones(5,1);

% % initial stage 1 generated from initial_transfer5.mat
% tf = 1.046833273939775e+04;
% x0_i = [1.495746722510590;0.000001002669660;0.006129720493607;0.000000302161724;-0.000899607989820;-0.000000013286327];
% xcf = [1.495746722510590;0.000001002669660;0.006129720493607;0.000000302161724;-0.000899607989820;-0.000000013286327];
% h0_i = 1e-3*ones(6,1);
% beta_i = 1e-3*ones(5,1);

x0_i = [1.9716164836809582;-0.0000016982488024;0.0109288701219422;0.0000961012751535;-0.0011435009582226;-0.0000207112628710];
h0_i = [-0.0013074412970034;-0.0005974031433551;0.0008808988541565;-0.1003654718445301;-0.0000002464425447;0.0000550503385134];
beta_i = [0.1000000000000000;-0.2843912220700550;-0.6394506321576822;-0.3239923337577765;0.1003282655938541];
tf = 12305.6202296647479670;
xcf = [3.8210482495396714;-0.0000000000000325;-0.2176791030142385;0.0001665010648302;-0.0017919029209414;-0.0000159320533486];


%% CALL THE SHOOTING METHOD

[sol_output] = asteroid_shooting_z(x0_i,h0_i,beta_i,xcf,tf,constants);

end_time = toc(start_time);
fprintf('FINISHED ASTEROID SHOOTING IN %5.4f sec\n',end_time);
