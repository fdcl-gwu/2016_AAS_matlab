% 7 Mar 16
% test gradient by comparing to numerical differentiation
clc
clear all
close all

% load constants and gravity model
constants = load_constants('castalia','false');
% perform preliminary computations for polyhedron potential model
constants.asteroid_grav = polyhedron_shape_input(constants);

tspan = [0,10000*1];
RelTol = constants.RelTol;
AbsTol = constants.AbsTol;

ode_options = odeset('RelTol',RelTol,'AbsTol',AbsTol);
constants.pot_model = 'mascon';

% pick a known position
pos = [1;1;1];
vel = rand(3,1);

state = [pos;vel];
% apply a small perturbation
deltapos = 1e-6*rand(size(pos));
deltavel = 1e-6*rand(size(vel));

deltastate = [deltapos;deltavel];

% calculate f(x + deltax) - f(x)
switch constants.asteroid_grav.num_f
    case 1024
        [U,dUdx, ~,~] = polyhedron_potential_mex_1024(pos, constants.asteroid_grav);
        [Udx,~, ~,~] = polyhedron_potential_mex_1024(pos+deltapos, constants.asteroid_grav);
    case 4092
        [U,dUdx, ~,~] = polyhedron_potential_mex_4092(pos, constants.asteroid_grav);
        [Udx,~, ~,~] = polyhedron_potential_mex_4092(pos+deltapos, constants.asteroid_grav);
end
% compare to gradient of f
diff_U = Udx - U;
Ugrad_diff = abs(dUdx'*deltapos - diff_U)

% test the mascon potential function
[U_mascon,dUdx_mascon] = mascon_potential(pos,constants.asteroid_grav,constants);
[Udx_mascon,~] = mascon_potential(pos+deltapos,constants.asteroid_grav,constants);

diff_Umascon = Udx_mascon-U_mascon;
Ugradd_diff_mascon = abs(dUdx_mascon'*deltapos - diff_Umascon)
% test the gradient of Jacobi integral
[J, dJdx] = jacobi(state,constants);
[Jdx, ~] = jacobi(state+deltastate,constants);

diff_J = Jdx - J;
Jgrad_diff = abs(dJdx'*deltastate - diff_J)

% test the jacobian of the state equations of motion (linearization)
f = ast_eoms(0,state,constants);
fdx = ast_eoms(0,state+deltastate,constants);

[dfdx, ~,~,~,~] = linearized_mat(state,constants.asteroid_grav,constants);

diff_f = fdx - f;
fgrad_diff = abs(dfdx*deltastate - diff_f)

%% plot trajectory and energy level during motion
% find equilibrium points
x_g=equilibria_finder(3,constants,constants.asteroid_grav);
% define initial state
% call ode45

% pick an initial condition
initial_state = [1.323575463621905;0;0.008105114223604;0.000000526662756;-0.000845286178738;-0.000000032580056];
[t,state]=ode113(@(t,state)ast_eoms(t,state,constants),tspan, initial_state,ode_options);

Jtraj = zeros(length(t),1);
for ii = 1:length(t)
    %calculate J over trajectory
    [Jtraj(ii), ~] = jacobi(state(ii,:)',constants);
end
traj_fig = figure;
hold all
grid on

vertex_plotter(constants.F,constants.V,traj_fig)
plot(state(:,1),state(:,2),'k')

figure
plot(t,Jtraj)
title('Jacobi')
ylabel('J')
xlabel('sec')

%% Test the STM computation over the trajectory
% define initial state
phi0 = eye(6,6);
phi0=reshape(phi0,36,1);
initial_stm = [initial_state;phi0];
%propogate the STM over the time span
[t,stm] = ode113(@(t,state)ast_stm(t,state,constants),tspan,initial_stm,ode_options);	% integration

% compute determinant at every time step
det_stm = zeros(1,length(t));
for ii =1:length(t)
    phi = reshape(stm(ii,7:end),6,6);
    det_stm(ii) = det(phi);
end

figure
plot(t,det_stm)
title('Det STM')
xlabel('sec')
ylabel('Det')