function [x0_out, cross_time_out,phi_cross_out, monodromy_matrix] = periodic_orbit_diffcorr(initial_state,constants)

TSPAN = constants.periodic_tspan ;

RelTol = constants.RelTol;
AbsTol = constants.AbsTol;
tol = constants.tol;
alpha = 0.01;
    
x0 = initial_state;

% logic to change the monodromy matrix based on the Poincare section
switch constants.periodic_diffcorr_section 
    case 'y_axis' % x = 0
        
    case 'x_axis' % y = 0
        rows_ind = [1 3 4 6];
        cols_ind = rows_ind;
        sec_ind = 2;
        secdot_ind = 5;
    case 'z_axis' % y = 0 or x = 0
        
end
% initial jacobi integral value
[J0, dJdx] = jacobi(x0,constants);

if constants.diffcorr_plot == 1
    fig_handle = figure;
    hold all;grid on
    title('Differential Correction Iterations')
    xlabel('X Axis')
    ylabel('Y Axis')
end

options_cross = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events',@(t,x)events_xcross(t,x,constants));
options_stm = odeset('RelTol',RelTol,'AbsTol',AbsTol);

delta = 1;
iter = 1;
while norm(delta) > tol && iter < 500
    fprintf('DiffCorr - %2g \n', iter)
    
    % propogate state and STM to first Poincare section crossing (do both
    % together)
    phi0 = eye(6,6);
    phi0=reshape(phi0,36,1);
    
    initial_stm = [x0;phi0];
    [t,state,cross_t,cross_state,ie] = ode113(@(t,state)ast_stm(t,state,constants),TSPAN,initial_stm,options_cross) ;

    % switch to plot successive iterations on the intial condition
    
    if constants.diffcorr_plot==1,
        set(0,'CurrentFigure',fig_handle)
        plot3(state(:,1),state(:,2),state(:,3));
        
        plot3(state(end,1),state(end,2),state(end,3),'bo');
        drawnow
    end
    
    % state at Poincare section crossing - logic to determine which of the
    % crossings to choose
    x1 = cross_state(end,1:6)';
    % calculate the state derivatives (f(x) at the crossing time)
    x1dot=ast_eoms(cross_t,x1,constants);
    
    % pull out the final stm (phi(t1,t0) )
    phi1 = cross_state(end,7:end);
    phi1 = reshape(phi1,6,6);
    
    [~, dJdx0] = jacobi(x0,constants);
    [~, dJdx1] = jacobi(x1,constants);
    
    % form the monodromy matrix
    phi1_reduced = phi1(rows_ind,cols_ind);
    x1dot_reduced = x1dot(rows_ind);
    phi1_sec_reduced = phi1(sec_ind,cols_ind);
    
    secdot0 = x0(secdot_ind);
    secdot1 = x1(secdot_ind);
    
    phi1_reduced_col = phi1(rows_ind,secdot_ind);
    phi1_sec_secdot = phi1(sec_ind,secdot_ind);

    dJdx0_reduced = dJdx0(cols_ind)';
    
    % maps updates from \delta [x0 z0 xd0 zd0] to \delta [x1 z1 xd1 zd1]
    % for y = 0 section
    monodromy_mat = phi1_reduced - x1dot_reduced/secdot1*phi1_sec_reduced - ...
        1/secdot0*(phi1_reduced_col-1/secdot1*x1dot_reduced*phi1_sec_secdot) * ...
        dJdx0_reduced;

    delta = (eye(4,4)-monodromy_mat)\(x1(rows_ind)-x0(rows_ind));
    
    % pick out the update to apply to the next initial state

    % vary x 
    del = [delta(1);0;delta(2);0;0;0];

    x0 = x0 + del;
    
    iter = iter+1;
    
    
end % end of while loop to correct the initial condition

x0_out = x0;
cross_time_out = cross_t(end);
phi_cross_out = phi1;
monodromy_matrix = monodromy_mat;

fprintf('DiffCorr Finished \n')
end



