
tspan = [0, 3.154e7];
h_ref = 86400; % 1 day (in seconds)
t_list = linspace(tspan(1),tspan(2),num_steps+1);


[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(gravity_rate_func_tabby,tspan,X0,h_ref,orbit_params)

