function [con, jac] = objective_z(xg, Q, x_0,phi_d,constants)
% output the constraint and jacobian
num_seg = constants.num_seg;
num_mid = num_seg - 1;
num_states = constants.num_states;
num_steps = constants.num_steps;
num_con = constants.num_con;

t = constants.t;
% xg is a big vector of all the unknowns to be solved for
% unpack xg
h0 = xg(1:num_states);
beta = xg(length(xg)-num_con+1:end);
% parse out xm, hm interior constraint points
xm = zeros(num_states,num_mid);
hm = zeros(num_states,num_mid);

for mid = 1:num_mid
    srow_idx = num_states+num_states*2*(mid-1)+1;
    erow_idx8 = srow_idx + 2*num_states - 1;
    
    mid_comb = xg(srow_idx:erow_idx8);
    xm(:,mid) = mid_comb(1:num_states);
    hm(:,mid) = mid_comb(num_states+1:end);
    
end

% propogate state, costate, sensitivities for each segment
% Generate a new complete trajectory

constants.control_switch = 'on';

% propogate each segment
% x_i, lambda_i
x_i = zeros(num_steps/num_seg,num_states,num_seg); % segment trajectories
h_i = zeros(num_steps/num_seg,num_states,num_seg);

for seg = 1:num_seg
    % logic to propagate each segment
    if seg == 1
        [~,z] = asteroid_optimal_ode(x_0,h0,0,t(1,end),num_steps/num_seg,constants);
    else
        [~,z] = asteroid_optimal_ode(xm(:,seg-1),hm(:,seg-1),t(seg,1),t(seg,end),num_steps/num_seg,constants);
    end
    
    % store the propogated state/costate/sensitivites into arrays for each
    % segment
    x_i(:,:,seg) = z(:,1:6);
    h_i(:,:,seg) = z(:,7:12);
    
end

% form constraint and jacobian by looping over the stages

con = zeros(num_con+num_states+num_mid*2*num_states,1);
jac = zeros(num_con+num_states+num_mid*2*num_states, num_con+num_states+num_mid*2*num_states);

% loop over the interior segments
for seg = 2:num_seg-1
    % indicies for forming the jacobian and constraint equations
    srow_idx = num_states*2*(seg-1)+1;
    erow_idx = srow_idx + 2*num_states - 1;
    
    scol_idx = num_states+ 2*num_states*(seg-2)+1;
    ecol_idx = scol_idx + 4*num_states - 1;
    
    [con_seg] = seq_newton_z(seg, x_i, h_i,Q, x_0,xm, hm,beta, h0, phi_d,constants); % compute stage con and jacobian
    % combine into a large g, dgdx vectors for output
    con(srow_idx:erow_idx) = con_seg;
end

% compute values for the first stage

[con_seg] = seq_newton_z(1, x_i, h_i,Q, x_0,xm, hm, beta, h0,phi_d,constants);
con(1:2*num_states) = con_seg;


% compute for last stage
[con_seg] = seq_newton_z(num_seg, x_i, h_i,Q, x_0,xm, hm, beta, h0,phi_d,constants);

con(length(con)-(num_states+num_con-1):end) = con_seg;

end % end of objective function