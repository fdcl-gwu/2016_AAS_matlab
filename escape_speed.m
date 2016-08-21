% 27 April 16
% escape speed from body of asteroid

% load asteroid
clc
close all

constants = load_constants('castalia','false'); % true = 1024
asteroid_grav = polyhedron_shape_input(constants);
constants.asteroid_grav = asteroid_grav;


% extract the normal vectors to each face
normal_face = asteroid_grav.normal_face;
center_face = asteroid_grav.center_face;

V1 = asteroid_grav.V1;
num_f = asteroid_grav.num_f;

omega = constants.omega*[0;0;1];

% escape speed for each face
radius_center = zeros(num_f,1);
long = zeros(num_f,1);
lat = zeros(num_f,1);

v_escape = zeros(num_f,1);
hwaitbar = waitbar(0,'Calculating...');
for ii = 1:num_f
    rc = center_face(ii,:)';
    nhat = normal_face(ii,:)';
    
    [Ur,~, ~,~] = polyhedron_potential_mex_4092(rc, asteroid_grav);
    %     [Ur,~, ~,~] = polyhedron_potential(rc, asteroid_grav);
    
    if Ur == 0
        rc = rc + 0.000001*rc/norm(rc);
        [Ur,~, ~,~] = polyhedron_potential_mex_4092(rc, asteroid_grav);
    end
    
    Umax = max(Ur,constants.mu/norm(rc));
    
    % calculate lat,long,r at current position
    
    radius_center(ii) = sqrt(rc(1).^2 + rc(2).^2 + rc(3).^2);
    long(ii) = atan2(rc(2),rc(1)) * 180/pi;
    lat(ii) = asin(rc(3)./radius_center(ii)) * 180/pi;
    
    % escape speed
    v_esc =-(nhat'*cross(omega,rc)) + sqrt((nhat'*cross(omega,rc))^2 -cross(omega,rc)'*cross(omega,rc) + 2*Umax);
    v_escape(ii) = v_esc;
    
    waitbar(ii/num_f);
end

close(hwaitbar);
%% create height map contour plot

% interpolate between the vertices
[longi,lati] = meshgrid(-180:0.5:180, -90:0.5:90);
vi = griddata(long,lat,v_escape*1e3,longi,lati);

% contour plot
figure
[c,h] = contourf(longi,lati,vi,[0.3:0.1:0.8]);
clabel(c,h);
colormap(jet)
colorbar
title('Escape Speed m/sec','interpreter','latex');
xlabel('Longitude','interpreter','latex')
ylabel('Latitude','interpreter','latex')