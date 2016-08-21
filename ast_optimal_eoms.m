function [state_dot] = ast_optimal_eoms(t,state,constants)

state = real(state);

% unpack the state (x is 6x1 and h is 6x1)
pos = state(1:3);
vel = state(4:6);
lam_pos = state(7:9);
lam_vel = state(10:12);

% test for imaginary numbers

try
    
    switch constants.pot_model
        case 'polyhedron'
            switch constants.asteroid_grav.num_f
                case 1024
                    [U,U_grad,U_grad_mat, Ulaplace] = polyhedron_potential_mex_1024(pos, constants.asteroid_grav);
                case 4092
                    [U,U_grad,U_grad_mat, Ulaplace] = polyhedron_potential_mex_4092(pos, constants.asteroid_grav);
            end
        case 'mascon'
            [U,U_grad] = mascon_potential(state,constants.asteroid_grav,constants);
        case 'matlab'
            [U,U_grad, U_grad_mat, Ulaplace] = polyhedron_potential(pos, constants.asteroid_grav);
    end
    
catch ME
    fprintf('Error in polyhedron potential mex\n')
    %     keyboard
end

switch constants.control_switch
    case 'off'
        u_control = zeros(3,1);
    case 'on'
        % express control vector in terms of costates
        u_control = -constants.um *lam_vel/norm(lam_vel);
    case 'sub'
        rand_vec = constants.sub_opt;
        u_control = constants.um*rand_vec/norm(rand_vec);
    case 'min_control'
        if norm(lam_vel) < constants.um
            u_control = -lam_vel;
        elseif norm(lam_vel) > constants.um
            u_control = -constants.um *lam_vel/norm(lam_vel);
        end
end

% state derivatives
pos_dot = vel;
vel_dot = U_grad + [2*constants.omega*vel(2);-2*constants.omega*vel(1);0] + constants.omega^2*[pos(1);pos(2);0] + u_control;

lam_pos_dot = - (U_grad_mat + constants.omega^2*diag([1 1 0]))*lam_vel;
lam_vel_dot = - lam_pos - [0 2*constants.omega 0;-2*constants.omega 0 0;0 0 0]*lam_vel;

state_dot = [pos_dot;vel_dot;lam_pos_dot;lam_vel_dot];
end