% 17 Feb 2016
% equations of motion about an asteroid

function [state_dot] = rot_asteroid(t,state,constants)

% unpack the state
pos = state(1:3);
vel = state(4:6);
% unpack constants

% calculate potential
[U0,U2, U_grad] = grav_potential(state,constants);
% state derivatives
A = [0 2 0;-2 0 0;0 0 0];
pos_dot = vel;
vel_dot = [pos(1);pos(2);0] + A*vel - U0*pos +U_grad;

state_dot = [pos_dot;vel_dot];

