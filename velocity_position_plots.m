function velocity_position_plots()
    

% function that will plot and solve the local error for different
    % values of h_ref, and plot them against each other. Also uses loglog
    % regression to solve for k and p values of error
    
    tspan =  [1,10];
    X0 = solution01(tspan(1));
    h_ref = 0.1;

    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];



    % set up params for sun and planet
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    % initial conditions
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    V0 = [x0; y0; dxdt0; dydt0];


     my_rate_func = @(t_in,V_in) gravity_rate_func_tabby(t_in,V_in,orbit_params);
     my_solution_func = @(t_in) compute_planetary_motion(t_in,V0,orbit_params);

     [t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration ...
         (my_rate_func,tspan,my_solution_func(tspan(1)),h_ref,DormandPrince, 5, .0001);
 
    X_solution = my_solution_func(t_list);


error_threshold = 1e-4; % Example threshold, adjust as needed


% Position vs Time
figure;
% plot(t_list, X_list(:,1), 'ro-', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', 2);
xlabel('Time');
ylabel('Position (x)');
title('Position vs Time');

% Velocity vs Time
figure;
% plot(t_list, X_list(:,3), 'bo-', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', 2);
xlabel('Time');
ylabel('Velocity (dx/dt)');
title('Velocity vs Time');

% Plot 2: Scatter plot of step size as a function of distance from the sun
r_list = sqrt(X_list(:,1).^2 + X_list(:,2).^2); % Calculate distance r = sqrt(x^2 + y^2)
h_list = diff(t_list); % Step sizes from time differences

figure;
scatter(r_list(1:end-1), h_list, 'filled');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Distance (r)');
ylabel('Step Size (h)');
title('Step Size as a Function of Distance from Sun');

% Plot options
figure;
subplot(3,1,1);
semilogx(r_list(1:end-1), h_list, 'o-');
xlabel('Distance (r)');
ylabel('Step Size (h)');
title('Step Size vs Distance (semilogx)');

subplot(3,1,2);
semilogy(r_list(1:end-1), h_list, 'o-');
xlabel('Distance (r)');
ylabel('Step Size (h)');
title('Step Size vs Distance (semilogy)');

subplot(3,1,3);
loglog(r_list(1:end-1), h_list, 'o-');
xlabel('Distance (r)');
ylabel('Step Size (h)');
title('Step Size vs Distance (loglog)');

end

