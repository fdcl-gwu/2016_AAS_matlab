% 25 July 16
% code to analyze all of the stages
clc
clearvars
close all
%% ASTEROID CONSTANTS
constants = load_constants('castalia','true'); % only 1024 faces
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

constants.pot_model = 'polyhedron'; % or mascon or matlab

constants.ode_options = odeset('RelTol',1e-9,'AbsTol',1e-9);
constants.optfsolve = optimoptions(@fsolve,'Display','off','TolFun',1e-5,'TolX',1e-9,...
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

%% Initialize PLOTS
% generate and save the handles to a bunch of figures
num_figs = 7;
fontsize = 18;
fontname = 'Times';

fig_title = {'Poincar\`e Section','Poincar\`e Section', 'Reach Distance','Reach Distance','Reach Distance','Trajectory', 'Control'};
fig_xlabel = {'$x$','$z$','$\phi_1$','$\phi_2$','$\phi_3$','$x$', '$t (sec)$'};
fig_ylabel = {'$\dot x$','$\dot z$','$d$','$d$','$d$','$y$','$u$ (mm/sec)'};

fig_handle = zeros(num_figs);

for ii = 1:num_figs
    fig_handle(ii) = figure();
    set(0, 'CurrentFigure', fig_handle(ii)) 
    hold all
    grid on
    title(fig_title(ii),'interpreter','latex','FontName',fontname,'FontSize',fontsize);
    xlabel(fig_xlabel(ii),'interpreter','latex','FontName',fontname,'FontSize',fontsize);
    ylabel(fig_ylabel(ii),'interpreter','latex','FontName',fontname,'FontSize',fontsize);
    set(gca,'FontName',fontname,'FontSize',fontsize);
end

C = {'magenta','cyan','green','blue','yellow','red'};
%% load the initial transfer and target
load('./results/initial_transfer4.mat');

time_total = [];
control_total = [];
state_total = [];

% loop over and load each stage of transfer
for stage = 1:4
    % compute the minimum of each stage
    fname = sprintf('./results/4constraint_transfer/hpc_stage%d_um1e-7.mat',stage);
    load(fname)
    
    [min_reach,reach_struct] = minimum_reach(sol_output,state_target(1,:), 'min_dist');
        
    % plot each stage in a seperate color
    reach_state = cat(1,reach_struct(:).reach_end);
    pmap = sol_output(1).constants.pmap;
    
    phi_d = cat(1,reach_struct(:).phi_d)*180/pi;
    dist = [reach_struct(:).dist];
    min_phi_d = reach_struct(min_reach.index).phi_d*180/pi;
    min_dist = min_reach.dist;
    
    % Poincare section x vs xdot
    set(0,'CurrentFigure',fig_handle(1));
    plot(reach_state(:,pmap(1)),reach_state(:,pmap(3)),'Color',C{stage},'Marker','x','LineStyle','none')
    plot(reach_state(min_reach.index,pmap(1)),reach_state(min_reach.index,pmap(3)),'ro','MarkerSize',10);
    text(reach_state(min_reach.index,pmap(1)),reach_state(min_reach.index,pmap(3)),sprintf('Stage %d',stage),'interpreter','Latex');
    
    % Poincare section z vs zdot
    set(0,'CurrentFigure',fig_handle(2));
    plot(reach_state(:,pmap(2)),reach_state(:,pmap(4)),'Color',C{stage},'Marker','x','LineStyle','none')
    plot(reach_state(min_reach.index,pmap(2)),reach_state(min_reach.index,pmap(4)),'ro','MarkerSize',10);
    text(reach_state(min_reach.index,pmap(2)),reach_state(min_reach.index,pmap(4)),sprintf('Stage %d',stage),'interpreter','Latex');
    
    % Distance vs Phi1
    set(0,'CurrentFigure',fig_handle(3));
    plot(phi_d(:,1),dist,'Color',C{stage},'Marker','x','LineStyle','none');
    plot(min_phi_d(1),min_dist,'ro','markersize',10)
    text(min_phi_d(1),min_dist,sprintf('Stage %d',stage),'interpreter','Latex');
    
    % Distance vs Phi2
    set(0,'CurrentFigure',fig_handle(4));
    plot(phi_d(:,2),dist,'Color',C{stage},'Marker','x','LineStyle','none');
    plot(min_phi_d(2),min_dist,'ro','markersize',10)
    text(min_phi_d(2),min_dist,sprintf('Stage %d',stage),'interpreter','Latex');
    
    % Distance vs Phi3
    set(0,'CurrentFigure',fig_handle(5));
    plot(phi_d(:,3),dist,'Color',C{stage},'Marker','x','LineStyle','none');
    plot(min_phi_d(3),min_dist,'ro','markersize',10)
    text(min_phi_d(3),min_dist,sprintf('Stage %d',stage),'interpreter','Latex');
    
    % Trajectory
    set(0,'CurrentFigure',fig_handle(6));
    plot3(min_reach.traj(:,1),min_reach.traj(:,2),min_reach.traj(:,3),'Color',C{stage});
    
    % combine the time, state, and control histories all together
    if stage > 1
        time_total = [time_total;min_reach.time+time_total(end)];
        control_total = [control_total;min_reach.control];
        state_total = [state_total;min_reach.traj];
    elseif stage == 1
        time_total = [time_total;min_reach.time];
        control_total = [control_total;min_reach.control];
        state_total = [state_total;min_reach.traj];
    end
    
    
end

% plot the initial and target states on each of the plots
% Poincare section x vs xdot
set(0,'CurrentFigure',fig_handle(1));
plot(state_initial(end,pmap(1)),state_initial(end,pmap(3)),'ks','markersize',10)
text(state_initial(end,pmap(1)),state_initial(end,pmap(3)),'Initial','interpreter','Latex','HorizontalAlignment','right','VerticalAlignment','bottom');

plot(state_target(end,pmap(1)),state_target(end,pmap(3)),'ks','markersize',10)
text(state_target(end,pmap(1)),state_target(end,pmap(3)),'Target','interpreter','Latex','HorizontalAlignment','left','VerticalAlignment','bottom');

% Poincare section z vs zdot
set(0,'CurrentFigure',fig_handle(2));
plot(state_initial(end,pmap(2)),state_initial(end,pmap(4)),'ks','markersize',10)
text(state_initial(end,pmap(2)),state_initial(end,pmap(4)),'Initial','interpreter','Latex','HorizontalAlignment','left','VerticalAlignment','bottom');
plot(state_target(end,pmap(2)),state_target(end,pmap(4)),'ks','markersize',10)
text(state_initial(end,pmap(2)),state_initial(end,pmap(4)),'Target','interpreter','Latex','HorizontalAlignment','right','VerticalAlignment','bottom');

% Trajectory
set(0,'CurrentFigure',fig_handle(6));
vertex_plotter(sol_output(1).constants.F,sol_output(1).constants.V,fig_handle(6));
plot3(state_initial(:,1),state_initial(:,2),state_initial(:,3),'k');
plot3(state_target(:,1),state_target(:,2),state_target(:,3),'k');
zlabel('$z$','interpreter','latex','FontName',fontname,'FontSize',fontsize);

%% COMPUTE FINAL TRAJECTORY TO TARGET

% propogate the orbit until it intersects the x axis again
constants.periodic_diffcorr_section = 'x_axis';
options_cross = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol,'Events',@(t,x)events_xcross(t,x,constants));
num_steps = constants.num_steps;
[t,state,cross_t,cross_state,ie] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,30000,constants.num_steps),min_reach.xf(1:6),options_cross);

xt = state_target(1,:);
x0 = min_reach.xf(1:6)'; % final stage
tf = cross_t(end,:);

% [sol_output] = asteroid_shooting_min(x0,xt,[0, tf],constants);
load('./results/4constraint_transfer/final_transfer.mat')

% plot the final leg
set(0,'CurrentFigure',fig_handle(6));
plot3(sol_output.traj(:,1),sol_output.traj(:,2),sol_output.traj(:,3),'r','linewidth',3)

time_total = [time_total;sol_output.t+time_total(end)];
control_total = [control_total;sol_output.control]*1e6;
state_total = [state_total;sol_output.state];

% compute the final control input
% Control
set(0,'CurrentFigure',fig_handle(7));
plot(time_total,control_total(:,1),'r')
plot(time_total,control_total(:,2),'g')
plot(time_total,control_total(:,3),'b')

h_leg = legend('$u_x$','$u_y$','$u_z$');
set(h_leg,'interpreter','Latex','FontName',fontname,'FontSize',fontsize);
