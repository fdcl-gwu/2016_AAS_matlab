% 10 March 2016

% find linearized EOMS for motion about an asteroid with a polyhedral
% potential field

function [A, U, Ug, Ugmat, Ulap] = linearized_mat(state,asteroid_grav,constants)

pos = state(1:3);
vel = state(4:6);

% calculate the potential at this state
switch constants.pot_model
    case 'polyhedron'
        switch constants.asteroid_grav.num_f
            case 1024
                [U,Ug, Ugmat, Ulap] = polyhedron_potential_mex_1024(pos, asteroid_grav);
            case 4092
                [U,Ug, Ugmat, Ulap] = polyhedron_potential_mex_4092(pos, asteroid_grav);
        end
    case 'mascon'
        [U,Ug,Ugmat] = mascon_potential(state,asteroid_grav,constants);
        Ulap = 0;
    case 'matlab'
        [U,Ug, Ugmat, Ulap] = polyhedron_potential(pos, constants.asteroid_grav);
end
A = zeros(6,6);

A(1:3,4:6) = eye(3,3);
A(4:6,1:3) = [constants.omega^2,0,0;0,constants.omega^2,0;0,0,0] + Ugmat;
A(4:6,4:6) = [0 2*constants.omega 0;-2*constants.omega 0 0;0 0 0];

