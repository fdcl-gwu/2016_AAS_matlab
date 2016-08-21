% 4 May 2016
% compare potential field models by running several tests
clc
close all

%% load asteroid data
constants = load_constants('castalia','false');
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

%% calculate equilibrium points
eqpoint_polyhedron = zeros(3,4);
eqpoint_mascon = zeros(3,4);
for ii = 1:4
    % eq point for polyhedron
    constants.pot_model = 'polyhedron';
    eqpoint_polyhedron(:,ii)=equilibria_finder(ii,constants,asteroid_grav);
    % eq point for mascon
    constants.pot_model = 'mascon';
    eqpoint_mascon(:,ii)=equilibria_finder(ii,constants,asteroid_grav);
end


%% contour plot of lines of constant gravitational attraction
% generate three figures for each axial contour plot
polyhedron_grav_fig_handles = zeros(1,3);
mascon_grav_fig_handles = zeros(1,3);
contour_planes = {'XY' 'XZ' 'YZ'};
for ii = 1:3
    polyhedron_grav_fig_handles(ii) = figure();
%     vertex_plotter(constants.F,constants.V,polyhedron_grav_fig_handles(ii))
    title(sprintf('Polyhedron Attraction %s',contour_planes{ii}),'interpreter','latex');
    
%     mascon_grav_fig_handles(ii) = figure();
%     %vertex_plotter(constants.F,constants.V,mascon_fig_handles(ii))
%     title(sprintf('Mascon Attraction %s',contour_planes{ii}),'interpreter','latex');
end

for ii = 1:3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % need to change the view of each plot for contours to work properly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    constants.plane = contour_planes{ii};
    
%     constants.pot_model = 'mascon';
%     grav_contour(mascon_grav_fig_handles(ii),constants);
    
    constants.pot_model = 'polyhedron';
    grav_contour(polyhedron_grav_fig_handles(ii),constants);
    
    
end

%% contour plot of potential
% generate three figures for each axial contour plot
polyhedron_pot_fig_handles = zeros(1,3);
mascon_pot_fig_handles = zeros(1,3);
contour_planes = {'XY' 'XZ' 'YZ'};
for ii = 1:3
    polyhedron_pot_fig_handles(ii) = figure();
    %     vertex_plotter(constants.F,constants.V,polyhedron_fig_handles(ii))
    title(sprintf('Polyhedron Potential %s',contour_planes{ii}),'interpreter','latex');
    
%     mascon_pot_fig_handles(ii) = figure();
    %     vertex_plotter(constants.F,constants.V,mascon_fig_handles(ii))
%     title(sprintf('Mascon Potential %s',contour_planes{ii}));
end

for ii = 1:3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % need to change the view of each plot for contours to work properly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    constants.plane = contour_planes{ii};
    
%     constants.pot_model = 'mascon';
%     pot_contour(mascon_pot_fig_handles(ii),constants);
    
    constants.pot_model = 'polyhedron';
    pot_contour(polyhedron_pot_fig_handles(ii),constants);

    
end

% %% contour plot of the difference between the two attraction models
% % generate three figures for each axial contour plot
% gravdiff_fig_handles = zeros(1,3);
% 
% contour_planes = {'xy' 'xz' 'yz'};
% for ii = 1:3
%     gravdiff_fig_handles(ii) = figure();
%     %     vertex_plotter(constants.F,constants.V,polyhedron_fig_handles(ii))
%     title(sprintf('Difference of Attraction %s',contour_planes{ii}));
% 
% end
% 
% for ii = 1:3
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % need to change the view of each plot for contours to work properly
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     constants.plane = contour_planes{ii};
%     
%     constants.pot_model = 'diff';
%     grav_contour(gravdiff_fig_handles(ii),constants);
% 
%     
%     
% end

% %% contour plot of difference of potential
% % generate three figures for each axial contour plot
% potdiff_fig_handles = zeros(1,3);
% contour_planes = {'xy' 'xz' 'yz'};
% for ii = 1:3
%     potdiff_fig_handles(ii) = figure();
%     %     vertex_plotter(constants.F,constants.V,polyhedron_fig_handles(ii))
%     title(sprintf('Difference of Potential %s',contour_planes{ii}));
% 
% end
% 
% for ii = 1:3
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % need to change the view of each plot for contours to work properly
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     constants.plane = contour_planes{ii};
%     
%     constants.pot_model = 'diff';
%     pot_contour(potdiff_fig_handles(ii),constants);
%     
%    
% end
% % zero velocity curves around the body
% 


% %% plot of direction of gradient of potential
% % generate three figures for each axial contour plot
% anglediff_fig_handles = zeros(1,3);
% contour_planes = {'xy' 'xz' 'yz'};
% for ii = 1:3
%     anglediff_fig_handles(ii) = figure();
%     %     vertex_plotter(constants.F,constants.V,polyhedron_fig_handles(ii))
%     title(sprintf('Angle of gradient difference %s',contour_planes{ii}));
% 
% end
% 
% for ii = 1:3
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % need to change the view of each plot for contours to work properly
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     constants.plane = contour_planes{ii};
%     
%     angle_contour(anglediff_fig_handles(ii),constants);
%     
%    
% end

%% Plot the zero velocity surfaces 
% generate three figures for each view
polyhedron_zvc_fig_handles = zeros(1,3);
mascon_zvc_fig_handles = zeros(1,3);
contour_planes = {'XY' 'XZ' 'YZ'};
for ii = 1:3
    polyhedron_zvc_fig_handles(ii) = figure();
    %     vertex_plotter(constants.F,constants.V,polyhedron_fig_handles(ii))
    title(sprintf('Polyhedron Zero-velocity Contours %s',contour_planes{ii}),'interpreter','latex');
    
%     mascon_zvc_fig_handles(ii) = figure();
    %     vertex_plotter(constants.F,constants.V,mascon_fig_handles(ii))
%     title(sprintf('Mascon Zero-velocity Contours %s',contour_planes{ii}),'interpreter','latex');
end

for ii = 1:3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % need to change the view of each plot for contours to work properly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constants.plane = contour_planes{ii};
    
%     constants.pot_model = 'mascon';
%     zvc_contour(mascon_grav_fig_handles(ii),constants);
    
    constants.pot_model = 'polyhedron';
    zvc_contour(polyhedron_zvc_fig_handles(ii),constants);
end
