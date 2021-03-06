% 26 July 2016
% ASTEROID TRANSFER DRIVER - now with extra z(tf) = 0 constraint
clc
close all
clearvars
%% ASTEROID CONSTANTS
constants = load_constants('castalia','true'); % only 1024 faces
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

constants.pot_model = 'polyhedron'; % or mascon or matlab

constants.ode_options = odeset('RelTol',1e-9,'AbsTol',1e-9);
constants.optfsolve = optimoptions(@fsolve,'Display','iter','TolFun',1e-5,'TolX',1e-9,...
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

constants.filename = 'second_stage_z.mat';
%% INITIALIZE TRANSFER
% this generates the initial and desired orbits (not any of the
% subintervals)
fprintf('INITIALIZE TRANSFER\n');

% constants.control_switch = 'off';
% [t_initial,state_initial,t_target,state_target]=initialize_transfer(constants);

load('./results/initial_transfer4.mat');
% 
% tf = t_initial(end);
% xcf = state_initial(end,1:6)';
% 
% % INITIAL GUESS OF COSTATE AND LAGRANGE MULTIPLIERS
% x0_i = state_initial(end,:)';
% h0_i = 1e-3*ones(6,1);
% beta_i = 1e-3*ones(5,1);

% fourth stage initial condition (from hpc_stage3_z.mat)
x0_i =[3.096001443589146;...
  -0.000002630129109;...
  -0.015716019220847;...
   0.000149692928676;...
  -0.001408643703713;...
   0.000141105337969];
h0_i =[0.000855757072297;...
  -0.001169268330287;...
  -0.000858465052047;...
  -0.100001782336235;...
   0.000003668063144;...
  -0.000126841493541];
beta_i =[0.100000000000000;...
   0.248031563826544;...
   0.178046006714276;...
  -0.080804108025452;...
   0.100020908809665];
tf = 1.410404717375811e+04;
xcf =[4.662931023131790;-0.000000000000239;1.881908402364113;0.000083719543268;-0.002053268068572;0.000125315421652];
%% CALL THE SHOOTING METHOD

% [sol_output] = asteroid_shooting_z(x0_i,h0_i,beta_i,xcf,tf,constants);
load('./results/5constraint_transfer/hpc_stage3_z_hightol.mat')
%% DETERMINE THE CLOSEST REACH STATE

[min_reach,reach_struct] = minimum_reach(sol_output,state_target(1,:), 'min_dist');

% do some optional plotting
plot_output(sol_output,reach_struct,min_reach);
%% INITIALIZE NEXT STAGE

% periodic orbit near the final state from the previous asteroid shooting
% constants.periodic_diffcorr_section = 'x_axis';
% constants.diffcorr_plot=1;
% [x0_out, cross_time_out, phi_cross_out,monodromy_matrix] = periodic_orbit_diffcorr(min_reach.xf(1:6)',constants);

% propogate the minimum orbit with no control and find the next positive x
% axis crossing

constants.periodic_diffcorr_section = 'x_axis';
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@(t,x)events_xcross(t,x,constants));
num_steps = constants.num_steps;
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,30000,constants.num_steps),min_reach.xf(1:6),options_cross);

% plot this orbit to visualize it
traj_fig = figure();
grid on
hold all
vertex_plotter(sol_output(1).constants.F,sol_output(1).constants.V,traj_fig);
plot3(state(:,1),state(:,2),state(:,3),'k')

% define the input data for the next stage
x0_i = min_reach.traj(end,:)';
h0_i = min_reach.costate(end,:)';
beta_i = sol_output(min_reach.index).beta;
tf = cross_t(end,:);
xcf = cross_state(end,:);
% call asteroid shooting again
% need to define x0_i, h0_i, beta_i, xcf,tf

% format and print the data to the command window
fprintf('x0_i = [%.16f;%.16f;%.16f;%.16f;%.16f;%.16f];\n',x0_i)
fprintf('h0_i = [%.16f;%.16f;%.16f;%.16f;%.16f;%.16f];\n',h0_i)
fprintf('beta_i = [%.16f;%.16f;%.16f;%.16f;%.16f];\n',beta_i)
fprintf('tf = %.16f;\n',tf)
fprintf('xcf = [%.16f;%.16f;%.16f;%.16f;%.16f;%.16f];\n',xcf)