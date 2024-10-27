function runner()
    
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

    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration_tabby ...
    (@rate_func01,tspan,X0,h_ref,DormandPrince)
    solution01(t_list);
    % X_list
    figure(1)
    clf
    hold on
    plot(t_list, X_list, "r")
    plot(t_list, solution01(t_list), "b")
    
end