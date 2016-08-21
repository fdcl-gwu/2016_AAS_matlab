%% driver for asteroid EOMS
% 16 Feb 16

clc
close all

constants = load_constants('castalia','true');
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;
constants.pot_model = 'polyhedron';

% find equilibrium points
x_g=equilibria_finder(3,constants,asteroid_grav);
% define initial state
% call ode45
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-6);

fontname = 'Times';
fontsize = 18;

%% Calculate periodic orbits
initial_state = [1.5;0;0;0;-0.0009;0]; % retrograde orbit initial condition
[t,state]=ode113(@(t,state)ast_eoms(t,state,constants),[0,20000], initial_state,ode_options);
% set(0,'CurrentFigure',traj_fig)
% plot(state(:,1),state(:,2),'b')
% plot ZVC for this jacobi level
[J1, ~] = jacobi(initial_state,constants);
% grav_contour(J,traj_fig,constants);

% periodic orbit at near this initial condition
constants.periodic_diffcorr_section = 'x_axis';
constants.diffcorr_plot=0;
[x0_out1, cross_time_out1, phi_cross_out1,monodromy_matrix1] = periodic_orbit_diffcorr(initial_state,constants);

initial_state = [3.82;0;0;0;-0.00148;0];
[t,state]=ode113(@(t,state)ast_eoms(t,state,constants),[0,20000], initial_state,ode_options);
% set(0,'CurrentFigure',traj_fig)
plot(state(:,1),state(:,2),'r')
% plot ZVC for this jacobi level
[J2, ~] = jacobi(initial_state,constants);
% grav_contour(J,traj_fig,constants);

% periodic orbit near this initial condition
constants.periodic_diffcorr_section = 'x_axis';
constants.diffcorr_plot=1;
[x0_out2, cross_time_out2, phi_cross_out2,monodromy_matrix2] = periodic_orbit_diffcorr(initial_state,constants);


%% Propogate both periodic orbits
numsteps = 1e3;
[t1,state1_body] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,4*cross_time_out1,numsteps),x0_out1,ode_options);
[t2,state2_body] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,4*cross_time_out2,numsteps),x0_out2,ode_options);

% save data to a mat file
t_initial = t1;
state_initial = state1_body;
t_target = t2;
state_target = state2_body;

% save('./results/initial_transfer5.mat','t_initial','state_initial','t_target','state_target')

% convert state to inertial frame and plot both
state1_inertial = zeros(size(state1_body));
state2_inertial = zeros(size(state2_body));

for ii = 1:length(t1)
   % calculate rotation matrix
   theta = constants.omega*t1(ii);
   
   Rb2i = ROT3(theta);
   
   pos1_inertial = Rb2i*state1_body(ii,1:3)';
   vel1_inertial = state1_body(ii,4:6)' + cross([0;0;constants.omega],state1_body(ii,1:3)');
   
   pos2_inertial = Rb2i*state2_body(ii,1:3)';
   vel2_inertial = state2_body(ii,4:6)' + cross([0;0;constants.omega],state2_body(ii,1:3)');
   
   state1_inertial(ii,:) = [pos1_inertial' vel1_inertial'];
   state2_inertial(ii,:) = [pos2_inertial' vel2_inertial'];
end

% setup figures for plotting
inertial_fig = figure;
hold all
grid on
vertex_plotter(constants.F,constants.V,inertial_fig)

body_fig = figure;
hold all
grid on
vertex_plotter(constants.F,constants.V,body_fig)

% inertial frame
set(0,'CurrentFigure',inertial_fig)
plot3(state1_inertial(:,1),state1_inertial(:,2),state1_inertial(:,3),'k')
plot3(state2_inertial(:,1),state2_inertial(:,2),state2_inertial(:,3),'k')

text(state1_inertial(1,1),state1_inertial(1,2),state1_inertial(1,3),'Initial','HorizontalAlignment','right','interpreter','latex')
text(state2_inertial(1,1),state2_inertial(1,2),state2_inertial(1,3),'Target','HorizontalAlignment','left','interpreter','latex')

% body frame
set(0,'CurrentFigure',body_fig)
plot3(state1_body(:,1),state1_body(:,2),state1_body(:,3),'k')
plot3(state2_body(:,1),state2_body(:,2),state2_body(:,3),'k')

text(state1_body(1,1),state1_body(1,2),state1_body(1,3),'Initial','HorizontalAlignment','right','interpreter','latex')
text(state2_body(1,1),state2_body(1,2),state2_body(1,3),'Target','HorizontalAlignment','left','interpreter','latex')

title('Transfer Trajectory','interpreter','latex','FontName',fontname,'FontSize',fontsize);
xlabel('$x$ (km)','interpreter','latex','FontName',fontname,'FontSize',fontsize);
ylabel('$y$ (km)','interpreter','latex','FontName',fontname,'FontSize',fontsize);
zlabel('$z$ (km)','interpreter','latex','FontName',fontname,'FontSize',fontsize);
set(gca,'FontName',fontname,'FontSize',fontsize);