% 7 March 2016
% locate equilibrium points
function x_g=equilibria_finder(eq_point,constants,asteroid_grav)

% % load constants and gravity model
% constants = load_constants('castalia');
% % perform preliminary computations for polyhedron potential model
% asteroid_grav = polyhedron_shape_input(constants);

% resonant radius
rs = constants.res_radius;

% initial guess for eq point
switch eq_point
    case 1
        x0 = [-0.04;-0.8;0]; % -C
    case 2
        x0 = [-0.05;0.7;0]; % +C
    case 3
        x0 = [0.9;0;0]; % +S
    case 4
        x0 = [-0.9;0;0]; % -S
end

% options = optimoptions(@fsolve,'Display','iter','Jacobian','on','DerivativeCheck','on','TolFun',1e-16,'TolX',1e-16);
% [pos,F,exitflag,output,JAC] = fsolve(@(pos)obj(pos,asteroid_grav,constants),x0,options);

% gradient descent
x_g = x0;
delta_x = 1;
tol = 1e-9;
alpha = 0.1;

while norm(delta_x) > tol
    [F,df] = obj(x_g,asteroid_grav,constants);
    delta_x = -inv(df)*F;
    
    x_g = x_g + alpha*delta_x;
    
end


end


function [F,df] = obj(pos,asteroid_grav,constants)
% objective function with jacobian
% switch mex function based on number of faces
switch constants.pot_model
    case 'polyhedron'
        switch constants.asteroid_grav.num_f
            case 1024
                [U,dUdx, ddUddx,~] = polyhedron_potential_mex_1024(pos, asteroid_grav);
            case 4092
                [U,dUdx, ddUddx,~] = polyhedron_potential_mex_4092(pos, asteroid_grav);
        end
    case 'mascon'
        [U,dUdx, ddUddx] = mascon_potential(pos, asteroid_grav,constants);
end
F = dUdx + [constants.omega^2*pos(1);constants.omega^2*pos(2);0];
F = F*1e6;
df = ddUddx + constants.omega^2*diag([1;1;0]);
df = df*1e6;
end