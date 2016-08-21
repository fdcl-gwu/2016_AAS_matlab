% 14 April 16 numerical continuation

function [x0_array, cross_time_array,phi1_array,monodromy_array] = jacobi_continuation(x0,period, phi1,monodromy,constants)

% inputs x0 - initial state of a periodic orbit

% input initial state and monodromy matrix of a periodic orbit
num_steps = 25;

% choose correct section pick ell
switch constants.periodic_diffcorr_section
    case 'x_axis' % y = 0
        rows = [1 3 4 6];
        cols = rows;
        
        sec = 2;
        secdot = 5;
    case 'y_axis' % x = 0
        rows = [2 3 5 6];
        cols = rows;
        
        sec = 1;
        secdot = 4;
    case 'z_axis' % x or y = 0
        rows = [1 2 4 5];
        cols = rows;
        
        sec = 3;
        secdot = 6;
end

x0_array = zeros(6,num_steps);
cross_time_array = zeros(num_steps);
phi1_array = zeros(6,6,num_steps);
monodromy_array = zeros(4,4,num_steps);

x0_array(:,1) = x0;
cross_time_array(1) = period;
phi1_array(:,:,1) = phi1;
monodromy_array(:,:,1) = monodromy;

delta_jacobi = -constants.jacobi_step;

% loop over and change the jacobi energy
for ii = 2:num_steps
    statedot=ast_eoms(cross_time_array(ii-1),x0_array(:,ii-1),constants);
    
    corr_mat = (eye(4,4)-monodromy_array(:,:,ii-1)) \ ( phi1_array(rows,secdot,ii-1) - statedot(rows)/statedot(sec)*phi1_array(sec,secdot,ii-1)) * (-1/2/x0_array(secdot,ii-1));
    delta = corr_mat*delta_jacobi;
    
    [J0, ~] = jacobi(x0_array(:,ii-1),constants);
    % compute new state with the new jacobi value (decrease in jacobi energy)
    
    tmp = x0_array(:,ii-1);
    tmp(cols) = tmp(cols)+delta;
    tmp(sec) = 0;
    
    velnew = jacobi_constraint(J0-constants.jacobi_step,tmp, constants);
    
    tmp(secdot) = sign(x0_array(secdot,ii-1))*velnew;
    
    % differential correction on this new initial condition
    constants.periodic_diffcorr_section = 'x_axis';
    [x0_array(:,ii), cross_time_array(ii), phi1_array(:,:,ii),monodromy_array(:,:,ii)] = periodic_orbit_diffcorr(tmp,constants);
end