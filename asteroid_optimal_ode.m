function [t,state] = asteroid_optimal_ode(x0,h0,t0,tf,num_steps,constants)

initial_condition = [x0;h0];

% replace this with a variational integrator at some point
[t,state] = ode113(@(t,state)ast_optimal_eoms(t,state,constants),linspace(t0,tf,num_steps),initial_condition,constants.ode_options);

end