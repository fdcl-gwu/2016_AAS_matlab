% gravitational potential about an asteroid 
% use some logic to determine model to implement
% 16 Feb 2016

function [U0,U2, U_grad] = grav_potential(state,constants)
% assume everything is already nondimensionalized

x = state(1);
y = state(2);
z = state(3);
% triaxial ellipsoid model
% load parameters from constants structure
C20 = constants.C20;
C22 = constants.C22;
mu = 1; % nondimensionalized

r = sqrt(x^2+y^2+z^2);

% calculate potential
U0 = mu / r;
U2 = -mu * C20 *(x^2+y^2-2*z^2)/(2*r^5) + 3*mu*C22*(x^2-y^2)/r^5;

% first order partial derivatives (EOMS)
dU2dx = -C20*x/r^5 + 5*C20*x/(2*r^7)*(x^2+y^2-2*z^2) + 6*C22*x/r^5 - 15*C22*x*(x^2-y^2)/r^7;
dU2dy = -C20*y/r^5 + 5*C20*y/(2*r^7)*(x^2+y^2-2*z^2) - 6*C22*y/r^5 - 15*C22*y*(x^2-y^2)/r^7;
dU2dz = 2 *C20*z/r^5 + 5*C20*z/(2*r^7)*(x^2+y^2-2*z^2) - 15*C22*z*(x^2-y^2)/r^7;
% second order derivatives (STM)

% third order derivatives (Variational Integrator/Optimal Costate EOMS)

U_grad = [dU2dx;dU2dy;dU2dz];


