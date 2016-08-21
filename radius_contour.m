% 26 April 16
% plot radius contours for asteroid

% load asteroid
clc
close all

constants = load_constants('castalia','false');
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;


%% create height map contour plot

% interpolate between the vertices
[longi,lati] = meshgrid(-180:0.5:180, -90:0.5:90);
ri = griddata(constants.long,constants.lat,constants.r,longi,lati);

% contour plot
figure
[c,h] = contour(longi,lati,ri,[0.3:0.1:0.8]);
clabel(c,h);
title('Radius Contour','interpreter','latex')
xlabel('Longitude (deg)','interpreter','latex')
ylabel('Latitude (deg)','interpreter','latex')