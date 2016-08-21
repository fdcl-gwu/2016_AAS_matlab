% 21 Mar 16
% find dependence of equilibrium point position and number of faces
clear all
clc
close all

% loop over different number of faces in polyhedral potential model
num_faces = [64,128,256,512,1024,2048,4092];

eq_points = zeros(3,4,length(num_faces));
run_time = zeros(2,4,length(num_faces));
for ii = 1:length(num_faces)
    fprintf('Eq points with %4d faces\n',num_faces(ii))
    % load constants and gravity model
    constants = load_constants('castalia');
    [constants.F,constants.V] = reducepatch(constants.F,constants.V,num_faces(ii));
    ef_start = tic;
    % perform preliminary computations for polyhedron potential model
    asteroid_grav = polyhedron_shape_input(constants);
    ef_end = toc(ef_start);
    % loop over the four eq positions and find the correct one
    for jj = 1:4
        eq_start = tic;
        x_g=equilibria_finder(jj,constants,asteroid_grav);
        eq_end = toc(eq_start);
        run_time(1,jj,ii)  =ef_end;
        run_time(2,jj,ii) = eq_end;
        
        eq_points(:,jj,ii) = x_g;
    end
    
    
end
% plot the eq points
figure

subplot(4,3,1)
semilogx(num_faces,squeeze(eq_points(1,1,:)),'ro-')
title('X Pos')
ylabel('-C')
subplot(4,3,4)
semilogx(num_faces,squeeze(eq_points(1,2,:)),'bo-')
ylabel('+C')
subplot(4,3,7)
semilogx(num_faces,squeeze(eq_points(1,3,:)),'go-')
ylabel('+S')
subplot(4,3,10)
semilogx(num_faces,squeeze(eq_points(1,4,:)),'co-')
ylabel('-S')
xlabel('Num Faces')

subplot(4,3,2)
semilogx(num_faces,squeeze(eq_points(2,1,:)),'ro-')
title('Y Pos')
ylabel('-C')
subplot(4,3,5)
semilogx(num_faces,squeeze(eq_points(2,2,:)),'bo-')
ylabel('+C')
subplot(4,3,8)
semilogx(num_faces,squeeze(eq_points(2,3,:)),'go-')
ylabel('+S')
subplot(4,3,11)
semilogx(num_faces,squeeze(eq_points(2,4,:)),'co-')
ylabel('-S')
xlabel('Num Faces')

subplot(4,3,3)
semilogx(num_faces,squeeze(eq_points(3,1,:)),'ro-')
title('Z Pos')
ylabel('-C')
subplot(4,3,6)
semilogx(num_faces,squeeze(eq_points(3,2,:)),'bo-')
ylabel('+C')
subplot(4,3,9)
semilogx(num_faces,squeeze(eq_points(3,3,:)),'go-')
ylabel('+S')
subplot(4,3,12)
semilogx(num_faces,squeeze(eq_points(3,4,:)),'co-')
ylabel('-S')
xlabel('Num Faces')
