% 18 July 2016
% find the closest state in sol_output to the desired target 

function [min_reach,reach_struct] = minimum_reach(sol_output,target_state, reach_switch)
constants = sol_output(1).constants;

%% REACHABILITY SET - set of states that approximate the final reachability set on the Poincare section
num_steps = sol_output(1).constants.num_steps;
num_seg = sol_output(1).constants.num_seg;
num_sims = constants.num_angles^3;
num_states = sol_output(1).constants.num_states;
pmap = constants.pmap;

um = constants.um;

% create the reachability set structure
reach_struct(num_sims) = struct('state',[],'costate',[],'reach_end',[], 'control',[],'time',[], 'phi_d',[], ...
    'dist',[],'dist_vert',[]);
for ii = 1:num_sims % loop over theta angles (poincare directions)
    state = zeros(num_steps,num_states);
    costate = zeros(num_steps,num_states);
    time = zeros(num_steps,1);
    % loop over the segments and combine trajectories into a big array
    for jj = 1:num_seg
        x_i = sol_output(ii).x_i;
        h_i = sol_output(ii).h_i;
        t_i = sol_output(ii).t;
        start_idx = (jj-1)*num_steps/num_seg+1;
        end_idx = start_idx-1+num_steps/num_seg;
        state(start_idx:end_idx,:) = x_i(:,:,jj);
        costate(start_idx:end_idx,:) = h_i(:,:,jj);
        time(start_idx:end_idx,1) = t_i(jj,:)';
    end
    
    
    reach_struct(ii).state= state;
    reach_struct(ii).costate = costate;
    reach_struct(ii).reach_end = [state(end,:) costate(end,:)];
    
    % difference btwn final state and target
    x1diff = target_state(pmap(1)) - reach_struct(ii).reach_end(pmap(1)); % x
    x2diff = target_state(pmap(2)) - reach_struct(ii).reach_end(pmap(2)); % z
    x3diff = target_state(pmap(3)) - reach_struct(ii).reach_end(pmap(3)); % xdot
    x4diff = target_state(pmap(4)) - reach_struct(ii).reach_end(pmap(4)); % zdot
    
    switch reach_switch
        case 'min_dist'
            reach_struct(ii).dist = sqrt(x1diff^2 + x2diff^2 + x3diff^2 + x4diff^2);
        case 'min_vert'
            reach_struct(ii).dist = sqrt(x1diff^2 + 10*x2diff^2 + x3diff^2 + x4diff^2);
        case 'min_z'
            reach_struct(ii).dist = sqrt(x1diff^2 + 1e9*x2diff^2 + x3diff^2 + 1e9*x4diff^2);
    end
    
    
    
    reach_struct(ii).control = um*costate(:,4:6)./repmat(sqrt(costate(:,4).^2+costate(:,5).^2 + costate(:,6).^2),1,3);
    reach_struct(ii).time = time;
    reach_struct(ii).phi_d = sol_output(ii).phi_d;
end

reach_poincare = cat(1,reach_struct(:).reach_end);

%% COMPARE ALL THE REACH STATES TO FIND THE MINIMUM TO THE TARGET

switch reach_switch
    case 'min_vert' % minimize the vertical components (z distance)
        % find the minimum and the row/column of the minimum distance
        [dist, x_ind] = min([reach_struct.dist]);
                
        % pull out the minimum reach state and corresponding minimum manifold
        xf = reach_poincare(x_ind,:);
        
        min_traj = [reach_struct(x_ind).state ];
        min_costate = reach_struct(x_ind).costate;
        min_control = reach_struct(x_ind).control;
        min_time = reach_struct(x_ind).time;
    case 'min_z'
        % find the minimum and the row/column of the minimum distance
        [dist, x_ind] = min([reach_struct.dist]);
                
        % pull out the minimum reach state and corresponding minimum manifold
        xf = reach_poincare(x_ind,:);
        
        min_traj = [reach_struct(x_ind).state ];
        min_costate = reach_struct(x_ind).costate;
        min_control = reach_struct(x_ind).control;
        min_time = reach_struct(x_ind).time;
    case 'min_dist'
        % find the minimum and the row/column of the minimum distance
        [dist, x_ind] = min([reach_struct.dist]);
                
        % pull out the minimum reach state and corresponding minimum manifold
        xf = reach_poincare(x_ind,:);
        
        min_traj = [reach_struct(x_ind).state ];
        min_costate = reach_struct(x_ind).costate;
        min_control = reach_struct(x_ind).control;
        min_time = reach_struct(x_ind).time;
    case 'min_x'
        
        % find the minimum (furthest to the left) final x position of the trajectory
        
        [~,x_ind] = min(reach_poincare(:,1));
        dist = norm(target_state(pmap)-reach_poincare(x_ind,pmap));
                
        
        xf = reach_poincare(x_ind,:);
        
        min_traj = [reach_struct(x_ind).state ];
        min_costate = reach_struct(x_ind).costate;
        min_control = reach_struct(x_ind).control;
        min_time = reach_struct(x_ind).time;
    case 'max_x'
        
        % find the maximum (furthest to the right) final x position of the trajectory
        
        [~,x_ind] = max(reach_poincare(:,1));
        dist = norm(target_state(pmap)-reach_poincare(x_ind,pmap));
        xf = reach_poincare(x_ind,:);
        
        min_traj = [reach_struct(x_ind).state ];
        min_costate = reach_struct(x_ind).costate;
        min_control = reach_struct(x_ind).control;
        min_time = reach_struct(x_ind).time;
        
end

min_reach.traj = min_traj;
min_reach.costate = min_costate;
min_reach.control = min_control;
min_reach.xf = xf;
min_reach.time = min_time;
min_reach.index = x_ind;
min_reach.dist = dist;

