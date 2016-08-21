function [g] = seq_newton(stage, x_i, h_i, Q, x_0,xm, hm, beta, h0,phi_d,constants)
num_seg = constants.num_seg;
num_mid = num_seg - 1;
num_states = constants.num_states;

% put in logic to determine the number of constraints
num_con = constants.num_con;

% logic if the stage is the first or last stage and to form g, dgdx
if stage == num_seg % last stage
    % calculate the constraint and jacobian of constraint
    xf = x_i(end,:,stage);
    hf = h_i(end,:,stage);
    
    % terminal constraints
    [m, dmdx] = constraints_nodenom(xf',phi_d,constants);
    % newton iteration
    % form the vector g
    g1 =  hf' + Q * ( xf' - constants.xcf) - dmdx' * beta;
    g2 = m;
    
    g = [g1;g2];
    
    
elseif stage == 1 % first stage
    
    % calculate the constraint and jacobian of constraint
    xf = x_i(end,:,stage);
    hf = h_i(end,:,stage);
    
    g1 = xf' - xm(:,stage);
    g2 = hf' - hm(:,stage);
    
    g = [g1;g2];
    
    
else % all the interior stages
    
    xminus = x_i(end,:,stage);
    hminus = h_i(end,:,stage);
    
    g1 = xminus' - xm(:,stage);
    g2 = hminus' - hm(:,stage);
    
    g = [g1;g2];
    
    
end

end % seq newton function