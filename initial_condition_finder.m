% 21 April 16
% iterativly try multiple initial conditions to see if there's a periodic
% orbit

clc
close all
clear all

constants = load_constants('castalia','true');
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;
constants.pot_model = 'polyhedron';

% %% retrograde equitorial orbit initial state
% initial_state = [1.4887;0;0;0;-0.0009;0];
% constants.periodic_diffcorr_section = 'x_axis';
% constants.diffcorr_plot=1;
% [x0_out, cross_time_out, phi_cross_out,monodromy_matrix] = periodic_orbit_diffcorr(initial_state,constants);
% [J, ~] = jacobi(initial_state,constants);
% 
% % normalize by angular velocity
% J_retronorm = J/constants.omega^2;

%% find equilibrium points
x_g=equilibria_finder(3,constants,asteroid_grav);
% define initial state
% call ode45
ode_options = odeset('RelTol',1e-9,'AbsTol',1e-9);

traj_fig = figure;
hold on
grid on
vertex_plotter(constants.F,constants.V,traj_fig)

x_i = 1.5;
C_i = 1.0833;

[Up,~,~,~] = polyhedron_potential_mex_1024([x_i,0,0]', constants.asteroid_grav);
[Um,~] = mascon_potential([x_i,0,0]',constants.asteroid_grav,constants);
norm_vel = sqrt(-2*(constants.omega^2*C_i-1/2*constants.omega^2*(x_i^2)-Up));

% interesting initial state
initial_state = [1.5;0;0;0;norm_vel;0];
% norm_vel = 0.000413092980545;
for ii = 0:5
    % pick an initial condition
    step = -0.000001;
    initial_state = [x_i;0;0;0;norm_vel+ii*step;0];
    % check the current jacobi constant
    [t,state]=ode45(@(t,state)ast_eoms(t,state,constants),[0,1.5e4], initial_state,ode_options);
    [J, ~] = jacobi(initial_state,constants);
    
    fprintf('J: %g\n', J)
    plot3(state(:,1),state(:,2),state(:,3))
    drawnow
end

