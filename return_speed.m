% 28 April 16
% guaranteed return speed
clc
close all
% load asteroid
constants = load_constants('castalia','false'); % true = 1024
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;

% find equilibrium points
eqpoint = 4; % or 4 for \pm S
constants.pot_model = 'polyhedron';
x_g=equilibria_finder(eqpoint,constants,asteroid_grav);
% calculate energy of \pm S equlibrium points
[J, ~] = jacobi([x_g;0;0;0],constants);

% plot zero velocity curves at this energy level
zvc_fig = figure;
hold all
grid on
vertex_plotter(constants.F,constants.V,zvc_fig)

% constants.plane = 'xy';
% grav_contour(J,zvc_fig,constants)

% calculate return velocity over entire body

% extract the normal vectors to each face
center_face = asteroid_grav.center_face;

num_f = asteroid_grav.num_f;

omega = constants.omega*[0;0;1];

% escape speed for each face
radius_center = zeros(num_f,1);
long = zeros(num_f,1);
lat = zeros(num_f,1);

v_return = zeros(num_f,1);
hwaitbar = waitbar(0,'Calculating...');
for ii = 1:num_f
    rc = center_face(ii,:)';
        
    [Ur,~, ~,~] = polyhedron_potential_mex_4092(rc, asteroid_grav);
    %     [Ur,~, ~,~] = polyhedron_potential(rc, asteroid_grav);
    
    if Ur == 0
        rc = rc + 0.000001*rc/norm(rc);
        [Ur,~, ~,~] = polyhedron_potential_mex_4092(rc, asteroid_grav);
    end
    
    V = 1/2*constants.omega^2*(rc(1)^2+rc(2)^2) + Ur;
    
    % calculate lat,long,r at current position
    
    radius_center(ii) = sqrt(rc(1).^2 + rc(2).^2 + rc(3).^2);
    long(ii) = atan2(rc(2),rc(1)) * 180/pi;
    lat(ii) = asin(rc(3)./radius_center(ii)) * 180/pi;
    
    % escape speed
    v_return(ii) = sqrt(2*(V-J));
        
    waitbar(ii/num_f);
end

close(hwaitbar);
%% create height map contour plot

% interpolate between the vertices
[longi,lati] = meshgrid(-180:0.5:180, -90:0.5:90);
vi = griddata(long,lat,real(v_return)*1e3,longi,lati);

% contour plot
figure
[c,h] = contourf(longi,lati,vi,[0:0.01:0.2]);
clabel(c,h);
colormap(jet)
colorbar
title('Return Speed m/sec','interpreter','latex');
xlabel('Longitude','interpreter','latex')
ylabel('Latitude','interpreter','latex')
