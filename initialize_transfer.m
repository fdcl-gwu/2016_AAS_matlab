function [t1,state1_body,t2,state2_body]=initialize_transfer(constants)
% generate the initial and final transfer orbits for the total transfer
ode_options = constants.ode_options;

%% Calculate periodic orbits
initial_state = [1.5;0;0;0;-0.0009;0]; % retrograde orbit initial condition

% periodic orbit at near this initial condition
constants.periodic_diffcorr_section = 'x_axis';
constants.diffcorr_plot=0;
[x0_out1, cross_time_out1, phi_cross_out1,monodromy_matrix1] = periodic_orbit_diffcorr(initial_state,constants);

initial_state = [1.8;0;0;0;-0.00098;0];

% periodic orbit near this initial condition
constants.periodic_diffcorr_section = 'x_axis';
constants.diffcorr_plot=0;
[x0_out2, cross_time_out2, phi_cross_out2,monodromy_matrix2] = periodic_orbit_diffcorr(initial_state,constants);

%% propogate the orbits
num_steps = constants.num_steps;
[t1,state1_body] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,cross_time_out1,constants.num_steps),x0_out1,ode_options);
[t2,state2_body] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,cross_time_out2,constants.num_steps),x0_out2,ode_options);


end