

function global_truncation_error()
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

    

    figure(1)
    subplot(2,1,1)
    hold on
    plot(t_list,X_list(:,1),'r','linewidth',2);
    plot(t_list,X_list(:,2),'b','linewidth',2);

    plot(t_list,X_solution(:,1),'k--','linewidth',2)
    plot(t_list,X_solution(:,2),'k--','linewidth',2)

    subplot(2,1,2)
    hold on
    plot(t_list,X_list(:,3),'r','linewidth',2);
    plot(t_list,X_list(:,4),'b','linewidth',2);

    plot(t_list,X_solution(:,3),'k--','linewidth',2)
    plot(t_list,X_solution(:,4),'k--','linewidth',2)

    [X_list(end,:)',X_solution(end,:)']

    desired_error_range = logspace(-8,-2,100);
    norm(X_list(end,:)'-X_solution(end,:)')

    h_list = zeros(length(desired_error_range),1);
    global_error_list = zeros(length(desired_error_range),1);

    for n = 1:length(desired_error_range)
        desired_error = desired_error_range(n);
        [t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration ...
        (my_rate_func,tspan,my_solution_func(tspan(1)),h_ref,DormandPrince, 5, desired_error);

        X_solution = my_solution_func(t_list);

        global_error_list(n) = norm(X_list(end,:)'-X_solution(end,:)');
        h_list(n) = h_avg;
    end

    figure(2);
    loglog(h_list,global_error_list,'ro','markerfacecolor','r','markersize',2)

    filter_params.min_xval = 1e-10;
    filter_params.max_xval = 1e-1;
    filter_params.min_yval = 1e-13;
    filter_params.max_yval = 1e-1;

    [p_global1, k_global1] = loglog_fit(h_list, global_error_list, filter_params);

    hold on
    loglog(h_list,k_global1*h_list.^p_global1,'b')

    p_global1
    k_global1

    % controlled variables for rate_func01

    % a vector of 100 equally spaced h_ref values from 10e-5 to 10e1 we
    % will solve our integration functions
%     num_h_vals = 100;
%      h_refs = logspace(-5, 1, num_h_vals);

% 
%     tf = tspan(2);
%     X0_solution = solution01(tf);
%     global_error = zeros(length(h_refs),1);
%    
%    
%     X_sol = zeros(length(h_refs),length(X0_solution));
%     analytical = zeros(length(h_refs), 1);
% 
%     desired_errors = logspace(0, 1, num_h_vals);
%     
% 
%     for i = 1:length(h_refs)
% 
%         [t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration ...
% (@rate_func01,tspan,X0,h_ref,DormandPrince, 4.9980, desired_errors(i));
% 
%     X_sol= X0_solution;
%     global_error(i) = norm(X_list(end,:)'- X_sol);
%  
% 
%     analytical(i) = norm(X_sol - X0);
% 
%     end
%         
%     % using provided log regression function to solve for each method's p
%     % and k values
%     
%     filter_params.min_xval = 1e-5;
%     filter_params.max_xval = 1e-1;
%     filter_params.min_yval = 1e-14;
%     filter_params.max_yval = 1e0;
%     
%     [p_global, k_global] = loglog_fit(h_refs, global_error, filter_params);
%     [p_analytical, k_analytical] = loglog_fit(h_refs, analytical, filter_params);
% 
%      [p_global_t, k_global_t] = loglog_fit(tspan, global_error, filter_params)
%     [p_analytical_t, k_analytical_t] = loglog_fit(tspan, analytical, filter_params)
% 
% 
%     % plotting errors on a log scale
%     clf;
% 
%     % plotting calculated data
%     
% %     loglog(h_refs, global_error, 'ro','markerfacecolor','r','markersize',2)
% %     hold on
% %     ylim([10^-15 10]);
% %     xlim([10^-6 10]);
% %     loglog(h_refs, analytical, 'go', 'markerfacecolor', 'g', 'markersize', 2)
% 
%     % and now solving lines of best fit
% %     fit_line_x = 10e-7:0.1:10e1;
% %     fit_line_global = k_global*fit_line_x.^p_global;
% %     fit_line_analytical = k_analytical*fit_line_x.^p_analytical;
% 
% 
%         fit_line_global_t = k_global_t*tspan.^p_global_t;
%     fit_line_analytical_t = k_analytical_t*tspan.^p_analytical_t;
%         loglog(tspan, fit_line_global_t, 'r--','markerfacecolor','r','markersize',1)
% 
%     % plotting lines of best fit
% %     loglog(fit_line_x, fit_line_global, 'r--','markerfacecolor','r','markersize',1)
% %     loglog(fit_line_x, fit_line_analytical, 'g-','markerfacecolor','r','markersize',1)
%     title("Global Truncation Error Plot: XB1, XB2 Vs. Analytical Methods and |XB1-XB2|")
%     xlabel("h values [-]")
%     ylabel("error [-]")
%     legend("Global error", "Analytical", "Global error", "Analytical")
       
    
end

%% rate_func01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end

