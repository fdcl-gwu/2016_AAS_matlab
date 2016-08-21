% 12 April 2016
% function to solve for a specific velocity given a desired jacobi integral
% value

function secdot = jacobi_constraint(J0,state_new, constants)

[Un,~, ~, ~] = polyhedron_potential(state_new(1:3), constants.asteroid_grav);

% switch based on which section we are looking at
switch constants.periodic_diffcorr_section
    case 'x_axis'
        secdot = sqrt(2*(J0 + 1/2*constants.omega^2*(state_new(1)^2 + state_new(2)^2) + Un - state_new(4)^2 - state_new(6)^2));
    case 'y_axis'
        secdot = sqrt(2*(J0 + 1/2*constants.omega^2*(state_new(1)^2 + state_new(2)^2) + Un - state_new(5)^2 - state_new(6)^2));
    case 'z_axis'
        secdot = sqrt(2*(J0 + 1/2*constants.omega^2*(state_new(1)^2 + state_new(2)^2) + Un - state_new(4)^2 - state_new(5)^2));
end
