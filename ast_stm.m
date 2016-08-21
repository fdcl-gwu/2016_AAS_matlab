% ODE function to propogate the state transition matrix about asteroids
% propogates both the nl equations and the state transition matrix
% 1 April 16
function state_dot = ast_stm(t,state,constants)

% unpack the state
pos = state(1:3);
vel = state(4:6);

phi = state(7:42);
phi = reshape(phi,6,6);

% jacobian of dynamics
[A, U, Ug, Ugmat, Ulap] = linearized_mat([pos;vel],constants.asteroid_grav,constants);

% state derivatives
pos_dot = vel;
vel_dot = Ug + [2 * constants.omega*vel(2);-2 * constants.omega*vel(1);0] + constants.omega^2*[pos(1);pos(2);0];

% STM
phi_dot = A*phi;
phi_dot = reshape(phi_dot,36,1);

state_dot = [pos_dot;vel_dot;phi_dot];
