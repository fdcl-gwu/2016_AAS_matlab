% 26 July 2016
% Multiple shooting to approximate the reachable set about an asteroid with
% an extra z constraint

function [sol_output]= asteroid_shooting_z(x0_i,h0_i,beta_i,xcf,tf,constants)

% MULTIPLE SHOOTING - DIVIDE UP STATE/COSTATE HISTORIES
num_steps = constants.num_steps;
num_seg = constants.num_seg; % should be cleanly divisble by the number of steps
t = zeros(num_seg, num_steps/num_seg);
x_i = zeros(num_steps/num_seg,6,num_seg); % time vs. dim vs. segment
h_i = zeros(num_steps/num_seg,6,num_seg);

constants.num_seg = num_seg;
constants.control_switch = 'off';
for ii = 1:num_seg
    t(ii,:) = linspace((ii-1)*tf/num_seg,ii*tf/num_seg,num_steps/num_seg);
    % propogate the no control trajectory to generate the interior guess points
    if ii == 1
        % propogate from x0, h0
        [~,state] = asteroid_optimal_ode(x0_i,h0_i,0,t(ii,end),num_steps/num_seg,constants);
        x_i(:,:,ii) = state(:,1:6);
        h_i(:,:,ii) = state(:,7:12);
    else
        [~,state] = asteroid_optimal_ode(x_i(end,:,ii-1),h_i(end,:,ii-1),t(ii-1,end),t(ii,end),num_steps/num_seg,constants);
        x_i(:,:,ii) = state(:,1:6);
        h_i(:,:,ii) = state(:,7:12);
    end
end
h = t(1,2)-t(1,1); % step size
constants.h = h;
constants.t = t;

num_states = size(x_i,2);
constants.num_states = num_states;

num_mid = num_seg-1;

% COST FUNCTION Q - MAPPING FOR POINCARE SECTION
Q = diag([1 0 1 1 0 1]);

% no control final state to maximize against
constants.xcf = xcf;

x0 = x0_i; % this is modified by the loop
% DEFINE INITIAL GUESSES OF COSTATE AND BETA
h01 =    [0.001025538587149;
    0.001029284562816;
    0.001019833570702;
    0.075688989125740;
    0.091705693377778;
    0.086407013656273];
beta1 = [0.100000000000000;
    0.10000000000000000000;
    0.100000000000000;
    0.100000000000000;
    0.100014058504762];

num_con = 5;
constants.num_con = num_con;

% all the interior patch points for the multiple shooting
xm = squeeze(x_i(end,:,1:end-1));
hm = squeeze(h_i(end,:,1:end-1));

xm0 = reshape(xm,num_states,num_seg-1);
hm0 = reshape(hm,num_states,num_seg-1);

%% LOOP OVER DIRECTION ON POINCARE SECTION
num_angles = constants.num_angles;

sim_space = [num_angles, num_angles, num_angles];
num_sims = prod(sim_space);

constants.control_switch = 'on';
sol_output(num_sims) = struct('phi_d',[],'xg_sol',[], 'fval',[],'exitflag',[],'output',[],...
    'x_i',[],'h_i',[],'xm',[],'hm',[],'h0',[],'beta',[],'diff_x',[],...
    'diff_h',[],'m',[],'costate_necc',[],'t',[],'constants',[],'x0',[],'time',[]);

fprintf('STARTING POINCARE LOOP\n');

% LOOP OVER ALL THE ANGLES
parfor sol_index = 1:num_sims
    % track time for each iteration
    start_time = tic;
    [phi1_ind, phi2_ind, phi3_ind] = ind2sub(sim_space,sol_index);
    
    phi_d = poincare_angles([phi1_ind, phi2_ind,phi3_ind],num_angles);
    
    % form the xg guess vector
    % [h0, xm, hm, beta]
    xg = zeros(num_states + num_con + num_mid*2*num_states,1);
    for mid = 1:num_mid
        srow_idx = num_states+num_states*2*(mid-1)+1;
        erow_idx8 = srow_idx + 2*num_states - 1;
        
        xg(srow_idx:erow_idx8) = [xm0(:,mid);hm0(:,mid)];
    end
    xg(1:num_states) = h01;
    xg(length(xg)-num_con+1:end) = beta1;
    
    % call fsolve
    [xg,fval,exitflag,output] = fsolve(@(xg)objective_z(xg,Q, x0, phi_d,constants),xg,constants.optfsolve);
    % gives solutions for the intermediary point - now propagate each for
    % the segment history
    [x_i, h_i, xm, hm, h0, beta] = prop_seg(xg,x0,constants);
    %plot_seg(t,x_i,h_i,xm,hm, x0, h0,constants);
    
    % calculate all the minus states
    xminus = squeeze(x_i(end,:,:));
    hminus = squeeze(h_i(end,:,:));
    % compare to the plus states
    diff_x = xminus(:,1:size(xm,2)) - xm;
    diff_h = hminus(:,1:size(hm,2)) - hm;
    
    % neccessary condition on terminal costate
    [m, dmdx] = constraints_nodenom_z(xminus(:,end),phi_d,constants);
    
    costate_necc = hminus(:,end) + Q * ( xminus(:,end) - constants.xcf) - dmdx' * beta;
    
    % save to output structure for later plotting
    sol_output(sol_index).phi_d = phi_d;
    sol_output(sol_index).xg_sol = xg;
    sol_output(sol_index).fval = fval;
    sol_output(sol_index).exitflag = exitflag;
    sol_output(sol_index).output = output;
    sol_output(sol_index).x_i = x_i;
    sol_output(sol_index).h_i = h_i;
    sol_output(sol_index).xm = xm;
    sol_output(sol_index).hm = hm;
    sol_output(sol_index).h0 = h0;
    sol_output(sol_index).beta = beta;
    sol_output(sol_index).diff_x = diff_x;
    sol_output(sol_index).diff_h = diff_h;
    sol_output(sol_index).m = m;
    sol_output(sol_index).costate_necc = costate_necc;
    sol_output(sol_index).t = t;
    sol_output(sol_index).x0 = x0;
    
    end_time = toc(start_time);
    sol_output(sol_index).time = end_time;
    
    fprintf('Iter: %7d in %7.3f sec\n',sol_index,end_time);

end

sol_output(1).constants = constants;

save(constants.filename,'sol_output')
end