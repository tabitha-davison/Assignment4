function planetary_plot_variable_step()

    % set up params for sun and planet
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    gravity_rate_func = @(t_in,V_in) gravity_rate_func_tabby(t_in, V_in, orbit_params);
    
    % initial conditions
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    V0 = [x0; y0; dxdt0; dydt0];
    
    % set up time range
    t_range = linspace(0, 30, 100);
    h_ref = 0.000001;
    
    % compute true planetary motion using ground truth function
    V_list_truth = compute_planetary_motion(t_range, V0, orbit_params);
    
    % define Butcher tableau for diff RK methods
    method = 'Bogacki';
    
    % bogacki butcher tableau
    Bogacki = struct();
    Bogacki.C = [0, 1/2, 3/4, 1];
    Bogacki.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
    Bogacki.A = [0, 0, 0, 0;
                 1/2, 0, 0, 0;
                 0, 3/4, 0, 0;
                 2/9, 1/3, 4/9, 0];

    p = 4;
    error_d = 0.001;

    % compute next step using explicit RK integration
    [~, V_list_rk, ~, ~] = ...
    explicit_RK_variable_step_integration(gravity_rate_func, ...
    [t_range(1), t_range(end)], V0, h_ref, Bogacki, p, error_d);

    V_list_rk = V_list_rk';


    % % ode45 integration with gravity_rate_func
    % [~, V_ode45] = ode45(@(t, V) my_rate_func(t, V), t_range, V0);

    % plot results
    figure;
    hold on;
    axis equal; 
    axis square;
    axis([-10, 15, -10, 15]);

    % plot sun
    plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, ...
        'DisplayName', 'Sun'); 

    % plot ground truth
    plot(V_list_truth(:, 1), V_list_truth(:, 2), 'k-', 'DisplayName', ...
        'Ground Truth');

    % plot results from RK methods
    plot(V_list_rk(:, 1), V_list_rk(:, 2), ...
        'g--', 'LineWidth', 1.5, 'DisplayName', method);

    % add labels
    legend('Location', 'Best');
    title('Planetary Motion: Ground Truth vs Bogacki RK Method');
    xlabel('X Position');
    ylabel('Y Position');
    hold off;

end