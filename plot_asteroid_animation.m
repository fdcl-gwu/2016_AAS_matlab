% 24 Aug 2016 - animate the rotation of an asteroid

close all
clearvars
clc

constants = load_constants('castalia','false'); % only 1024 faces
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

type = 'none';
animation_fname = '.\animation\castalia';

time = linspace(0,86400,1e4);

inertial_fig = figure('color','white','position',[100 100 1124 868]);
hold on
axis off
axis equal
ast_inertial = patch('Faces',constants.F,'Vertices',constants.V);
set(ast_inertial,'FaceLighting','gouraud','AmbientStrength',0.3,...
    'DiffuseStrength',0.8,'SpecularStrength',0.9,'BackFaceLighting','unlit','Facealpha',1,'Facecolor',[0.8 0.8 0.8], 'EdgeAlpha',0.0);
axis([-0.9 0.9 -0.9 0.9 -0.4 0.4]);
lightangle(90,0)
view(3)

switch type
    case 'gif'
        
        f = getframe;
        [im,map] = rgb2ind(f.cdata,256,'nodither');
    case 'movie'
        nFrames = length(time);
        vidObj = VideoWriter(strcat(animation_fname,'.avi'));
        vidObj.Quality = 100;
        vidObj.FrameRate = 30;
        open(vidObj);
    case 'images'
    
    otherwise
        
end


h_waitbar = waitbar(0,'Starting Animation...');
% plot both in another loop and save a video to create the animation
for ii = 1:10:length(time)
    % calculate rotation matrix
    theta = constants.omega*time(ii);
    
    Rb2i = ROT3(theta);

    % inertial frame
    set(0,'CurrentFigure',inertial_fig);
    nV = constants.V*Rb2i';
    set(ast_inertial,'Vertices',nV);
    drawnow
  
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
    
    waitbar(ii / length(time))
end

close(h_waitbar);

% Output the movie as an avi file
switch type
    case 'gif'
        
        
    case 'movie'
        close(vidObj);
        
    otherwise
        
end
