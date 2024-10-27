function local_truncation_error()
    % function that will plot and solve the local error for different
    % values of h_ref, and plot them against each other. Also uses loglog
    % regression to solve for k and p values of error
    
    tspan =  [0,10];
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

    % controlled variables for rate_func01
    t1 = tspan(1);
    X01_sol = solution01(t1);
    
    X0_solution = X01_sol;

    % a vector of 100 equally spaced h_ref values from 10e-5 to 10e1 we
    % will solve our integration functions
    num_h_vals = 100;
    h_refs = logspace(-5, 1, num_h_vals);

    % establishing the size of the vectors we will store the error for a
    % given h value
    XB1_errors = zeros(length(h_refs),1);
    XB2_errors = zeros(length(h_refs),1);
    XB1_XB2_errors = zeros(length(h_refs),1);
    
    X_sol = zeros(length(X0_solution), length(h_refs));
    analytical = zeros(length(h_refs), 1);

    for i = 1:length(h_refs)
        % X_next with the midpoint and euer method
        [XB1, XB2, ~] = RK_step_embedded(@rate_func01,t1,X0,h_ref,DormandPrince)
        
        % solving for the analytical solution
        X_sol(:, i) = solution01(t1 + h_refs(i));
        
        % solving for the error of the three methods
        XB1_errors(i) = norm(XB1- X_sol(:, i));
        XB2_errors(i) = norm(XB2- X_sol(:, i));
        XB1_XB2_errors(i) = norm(XB1-XB2);
        
        analytical(i) = norm(X_sol(:, i) - X0_solution);
    end
        
    % using provided log regression function to solve for each method's p
    % and k values
    
    filter_params.min_xval = 10e-4;
    filter_params.max_xval = 10e-2;
    
    [p_XB1, k_XB1] = loglog_fit(h_refs, XB1_errors, filter_params);
    [p_XB2, k_XB2] = loglog_fit(h_refs, XB2_errors, filter_params);
    [p_XB1_XB2, k_XB1_XB2] = loglog_fit(h_refs, XB1_XB2_errors, filter_params);
    [p_analytical, k_analytical] = loglog_fit(h_refs, analytical, filter_params);

    % plotting errors on a log scale
    clf;

    % plotting calculated data
    loglog(h_refs, XB1_errors, 'ro','markerfacecolor','r','markersize',2)
    hold on
    loglog(h_refs, XB2_errors, 'bo','markerfacecolor','r','markersize',2)
    loglog(h_refs, XB1_XB2_errors, 'ko','markerfacecolor','r','markersize',2)
    loglog(h_refs, analytical, 'go', 'markerfacecolor', 'g', 'markersize', 2)

        % % plotting lines of best fit
        % loglog(fit_line_x, fit_line_exp_midpoint, 'c-','markerfacecolor','r','markersize',1)
        % loglog(fit_line_x, fit_line_exp_euler, 'm-','markerfacecolor','r','markersize',1)
        % loglog(fit_line_x, fit_line_analytical, 'y-','markerfacecolor','r','markersize',1)
        % title("Local Truncation Error Plot: Midpoint, Forward Euler Vs. Analytical Methods")
        % xlabel("h values [-]")
        % ylabel("error [-]")
        % legend("Midpoint Data", "Euler Data", "Analytical Data", "Midpoint Fit",...
        %     "Euler Fit", "Analytical Fit", 'location', 'nw')
    
end