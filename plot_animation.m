% 17 Aug 2016 - function to create an animation of the complete transfer
% trajectory

close all
clearvars
clc

constants = load_constants('castalia','true'); % only 1024 faces
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

fig_flag = 'body';
type = 'movie';
animation_fname = '.\animation\body';


switch fig_flag
    case 'body'
        % initialize a figure for the body fixed and inertial transfer animations
        body_fig = figure('color','white','position',[100 100 1124 868]);
        axis off
        axis equal
        hold on
        ast_body = patch('Faces',constants.F,'Vertices',constants.V);
        set(ast_body,'FaceLighting','gouraud','AmbientStrength',0.5,...
            'Facealpha',1,'Facecolor','green', 'EdgeAlpha',0.5);
        axis([-7 7 -7 7 -2 2]);
        view(3)
    case 'inertial'
        
        inertial_fig = figure('color','white','position',[100 100 1124 868]);
        hold on
        axis off
        axis equal
        ast_inertial = patch('Faces',constants.F,'Vertices',constants.V);
        set(ast_inertial,'FaceLighting','gouraud','AmbientStrength',0.5,...
            'Facealpha',1,'Facecolor','green', 'EdgeAlpha',0.5);
        axis([-7 7 -7 7 -2 2]);
        view(3)
end

% load the transfers and save to a big array
time_total = [];
control_body = [];
state_body = [];
seg_color = {}; % color of the segment based on which portion of the trajectory it is a part of

C = {'magenta','cyan','green','blue','yellow','red'};

% initial orbit
load('./results/initial_transfer4.mat');

% propogate initial and final orbits for a long period of time
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-6);
constants.pot_model = 'polyhedron';
numsteps = 5000;
[t_initial,state_initial] = ode113(@(t,state)bw_ast_eoms(t,state,constants),-linspace(4*t_initial(end),0,numsteps),state_initial(end,:),ode_options);
t_initial = -flipud(t_initial);
state_initial = flipud(state_initial);
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


switch type
    case 'gif'
        
        f = getframe;
        [im,map] = rgb2ind(f.cdata,256,'nodither');
    case 'movie'
        nFrames = length(time_total);
        vidObj = VideoWriter(strcat(animation_fname,'.avi'));
        vidObj.Quality = 100;
        vidObj.FrameRate = 30;
        open(vidObj);
    case 'images'
    
    otherwise
        
end


h_waitbar = waitbar(0,'Starting Animation...');
% plot both in another loop and save a video to create the animation
for ii = 1:10:length(time_total)
    % calculate rotation matrix
    theta = constants.omega*time_total(ii);
    
    Rb2i = ROT3(theta);
    
    pos_inertial = Rb2i*state_body(ii,1:3)';
    vel_inertial = state_body(ii,4:6)' + cross([0;0;constants.omega],state_body(ii,1:3)');
    
    state_inertial(ii,:) = [pos_inertial' vel_inertial'];
    control_inertial(ii,:) = (Rb2i*control_body(ii,:)')';
    
    switch fig_flag
        case 'body'
            % body fixed frame
            set(0,'CurrentFigure',body_fig);
            plot3(state_body(ii,1),state_body(ii,2),state_body(ii,3),'Color',seg_color{ii},'Marker','.','Markersize',9, 'LineStyle','-')
            drawnow
        case 'inertial'
            % inertial frame
            set(0,'CurrentFigure',inertial_fig);
            nV = constants.V*Rb2i';
            set(ast_inertial,'Vertices',nV);
            plot3(pos_inertial(1),pos_inertial(2),pos_inertial(3),'Color',seg_color{ii},'Marker','.','Markersize',9,'LineStyle','-')
            drawnow
    end
    
    switch type
        case 'gif'
            
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            outfile = strcat(animation_fname,'.gif');
            
            % On the first loop, create the file. In subsequent loops, append.
            if ii==1
                imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
            else
                imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
            end
        case 'movie'
            writeVideo(vidObj,getframe(gca));
        case 'images'
            outfile =  sprintf('%s_%d.jpg',animation_fname,ii);
            frame = getframe(1);
            imwrite(frame2im(frame),outfile);
        otherwise
            
    end
    
    waitbar(ii / length(time_total))
end

close(h_waitbar);

% Output the movie as an avi file
switch type
    case 'gif'
        
        
    case 'movie'
        close(vidObj);
        
    otherwise
        
end
