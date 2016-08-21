% 17 Aug 2016 - function to create an animation of the complete transfer
% trajectory

close all
clearvars
clc

constants = load_constants('castalia','true'); % only 1024 faces
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

% initialize a figure for the body fixed and inertial transfer animations
body_fig = figure();
axis equal
hold on

ast_body = patch('Faces',constants.F,'Vertices',constants.V);

set(ast_body,'FaceLighting','gouraud','AmbientStrength',0.5,...
    'Facealpha',1,'Facecolor','green', 'EdgeAlpha',0.5);

inertial_fig = figure();
hold on
axis equal
ast_inertial = patch('Faces',constants.F,'Vertices',constants.V);
set(ast_inertial,'FaceLighting','gouraud','AmbientStrength',0.5,...
    'Facealpha',1,'Facecolor','green', 'EdgeAlpha',0.5);

% load the transfers and save to a big array
time_total = [];
control_body = [];
state_body = [];
seg_color = {}; % color of the segment based on which portion of teh trajectory it is a part of

C = {'magenta','cyan','green','blue','yellow','red'};

% initial orbit
load('./results/initial_transfer4.mat');

% propogate initial and final orbits for a long period of time
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-6);
constants.pot_model = 'polyhedron';
numsteps = 5000;
[t_initial,state_initial] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,4*t_initial(end),numsteps),state_initial(1,:),ode_options);
[t_target,state_target] = ode113(@(t,state)ast_eoms(t,state,constants),linspace(0,22*t_target(end),numsteps),state_target(1,:),ode_options);


time_total = [time_total;t_initial];
state_body = [state_body;state_initial];
control_body = [control_body;zeros(length(t_initial),3)];
seg_color = cell(length(t_initial),1);
seg_color(:) = {'k'};

% loop over and load each stage of transfer
for stage = 1:4
    % compute the minimum of each stage
    fname = sprintf('./results/4constraint_transfer/hpc_stage%d_um1e-7.mat',stage);
    load(fname)
    
    [min_reach,reach_struct] = minimum_reach(sol_output,state_target(1,:), 'min_dist');
    
    % combine the time, state, and control histories all together
    
    time_total = [time_total;min_reach.time+time_total(end)];
    control_body = [control_body;min_reach.control];
    state_body = [state_body;min_reach.traj];
    for ii = 1:length(min_reach.time)
       seg_color{end+1} = C{stage}; 
    end
    
    
end

load('./results/4constraint_transfer/final_transfer.mat')

time_total = [time_total;sol_output.t+time_total(end)];
control_body = [control_body;sol_output.control]*1e6;
state_body = [state_body;sol_output.state];
for ii = 1:length(sol_output.t)
   seg_color{end+1} = 'red' ;
end

% target orbit
time_total = [time_total;t_target+time_total(end)];
state_body = [state_body;state_target];
control_body = [control_body;zeros(length(t_target),3)];
for ii = 1:length(t_target)
   seg_color{end+1} = 'k'; 
end
% convert state to inertial frame and plot both
state_inertial = zeros(size(state_body));
control_inertial = zeros(size(control_body));

nFrames = length(time_total);
vidObj = VideoWriter('animation.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 8;
open(vidObj);
    
% plot both in another loop and save a video to create the animation
for ii = 1:50:length(time_total)
    % calculate rotation matrix
    theta = constants.omega*time_total(ii);
    
    Rb2i = ROT3(theta);
    
    pos_inertial = Rb2i*state_body(ii,1:3)';
    vel_inertial = state_body(ii,4:6)' + cross([0;0;constants.omega],state_body(ii,1:3)');
    
    state_inertial(ii,:) = [pos_inertial' vel_inertial'];
    control_inertial(ii,:) = (Rb2i*control_body(ii,:)')';
    
    % body fixed frame
    set(0,'CurrentFigure',body_fig);
    plot3(state_body(ii,1),state_body(ii,2),state_body(ii,3),'Color',seg_color{ii},'Marker','.','Markersize',10, 'LineStyle','-')
    drawnow
    
    % inertial frame
    set(0,'CurrentFigure',inertial_fig);
    nV = constants.V*Rb2i';
    set(ast_inertial,'Vertices',nV);
    plot3(pos_inertial(1),pos_inertial(2),pos_inertial(3),'Color',seg_color{ii},'Marker','.','Markersize',10,'LineStyle','-')
    drawnow
    
    writeVideo(vidObj,getframe(gca));
    
    
end

close(vidObj);
