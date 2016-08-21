function [x_i, h_i, xm, hm, h0, beta] = prop_seg(xg,x0,constants)
num_seg = constants.num_seg;
num_mid = num_seg - 1;
num_states = constants.num_states;
num_steps = constants.num_steps;
num_con = constants.num_con;

t = constants.t;
% xg is a big vector of all the unknowns to be solved for
% unpack xg
eye8x8 = eye(8,8);

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