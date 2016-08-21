% 26 July 2016
% Multiple shooting to minimize to the final target state

function [sol_output]= asteroid_shooting_min(x0,xt,tspan,constants)

% MULTIPLE SHOOTING - DIVIDE UP STATE/COSTATE HISTORIES
num_steps = constants.num_steps;
num_seg = constants.num_seg; % should be cleanly divisble by the number of steps

num_states = size(x0,2);
constants.num_states = num_states;

num_mid = num_seg-1;

% target state to go towards
constants.xt = xt;

% DEFINE INITIAL GUESSES OF COSTATE
% h0 =[0.000918240948946;...
%    0.016223240681880;...
%   -0.003506901650284;...
%   -7.802190564281815;...
%    7.031321476164608;...
%    1.347079783658622];

h0 = 1*ones(6,1);

constants.control_switch = 'on';
sol_output = struct('xg_sol',[], 'fval',[],'exitflag',[],...
    'x_i',[],'h_i',[],'xm',[],'hm',[],'h0',[],'diff_x',[],...
    'diff_h',[],'costate_necc',[],'t',[],'constants',[],'x0',[],'time',[]);
optfsolve = optimoptions(@fsolve,'Display','iter','TolFun',1e-5,'TolX',1e-9,...
    'MaxIter',5000,'MaxFunEvals',5000, 'Algorithm', 'trust-region-reflective','Jacobian','off',...
    'DerivativeCheck','off', 'FunValCheck','on');

% call fsolve
[xg,fval,exitflag,~] = fsolve(@(xg)obj(xg,x0,tspan,constants),h0,optfsolve);

% propogate the final trajectory
[t,state] = asteroid_optimal_ode(x0,xg,0,tspan(end),num_steps,constants);

sol_output.xg_sol = xg;
sol_output.fval = fval;
sol_output.exitflag = exitflag;
sol_output.h0 = xg;

sol_output.costate = state(:,7:12);
sol_output.state = state(:,1:6);

sol_output.control = constants.um*sol_output.costate(:,4:6)./repmat(sqrt(sol_output.costate(:,4).^2+sol_output.costate(:,5).^2 + sol_output.costate(:,6).^2),1,3);

sol_output.t = t;
sol_output.traj = state;
sol_output.x0 = x0;

end

function F = obj(h0,x0,tspan,constants)


% integrate the equations of motion for the given h_0
constants.control_switch = 'min_control';
[~,state] = asteroid_optimal_ode(x0,h0,0,tspan(end),constants.num_steps,constants);

% subtract the final state from the desired final state
F = ( (state(end,1:6)'- constants.xt') );


end

function [x_i, h_i, xm, hm, h0] = prop_seg_min(x0,h0,tconstants)
num_seg = constants.num_seg;
num_mid = num_seg - 1;
num_states = constants.num_states;
num_steps = constants.num_steps;
num_con = constants.num_con;


% propogate state, costate, sensitivities for each segment
% Generate a new complete trajectory

constants.control_switch = 'on';

% propogate each segment/stage using the linearized system to generate
% x_i, lambda_i
x_i = zeros(num_steps/num_seg,num_states,num_seg);
h_i = zeros(num_steps/num_seg,num_states,num_seg);
for seg = 1:num_seg
    % loop over each STM to calculate the state
    state = zeros(size(x_i(:,:,seg)));
    costate = zeros(size(h_i(:,:,seg)));
    if seg == 1
        [~,z] = asteroid_optimal_ode(x0,h0,0,t(1,end),num_steps/num_seg,constants);
        
    else
        [~,z] = asteroid_optimal_ode(xm(:,seg-1),hm(:,seg-1),t(seg,1),t(seg,end),num_steps/num_seg,constants);
        
    end
    
    
    x_i(:,:,seg) = z(:,1:6);
    h_i(:,:,seg) = z(:,7:12);
    
    
end

end
