% 12 May 16
% follows convention used by Scheeres
function [J, dJdx] = jacobi(state,constants)

pos = state(1:3);
vel = state(4:6);

% calculate potential at this state
switch constants.pot_model
    case 'polyhedron'
switch constants.asteroid_grav.num_f
    case 1024
        [U,dUdp,~,~] = polyhedron_potential_mex_1024(pos, constants.asteroid_grav);
    case 4092
        [U,dUdp,~,~] = polyhedron_potential_mex_4092(pos, constants.asteroid_grav);
end
    case 'mascon'
        [U,dUdp] = mascon_potential(state,constants.asteroid_grav,constants);
    case 'matlab'
        [U,dUdp, ~, ~] = polyhedron_potential(pos, constants.asteroid_grav);
end
% integral value
J = 1/2*constants.omega^2*(pos(1)^2+pos(2)^2) + U - 1/2*(vel'*vel) ;

% gradient of the jacobi integral

dJdx = [constants.omega^2*pos(1) + dUdp(1);...
        constants.omega^2*pos(2) + dUdp(2);...
        dUdp(3);...
        -vel(1);...
        -vel(2);...
        -vel(3)];
