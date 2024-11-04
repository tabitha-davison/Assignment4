function adaptive_step_ex()
    % Set up parameters for sun and planet
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    % Define gravity rate function
    gravity_rate_func = @(t_in,V_in) gravity_rate_func_tabby(t_in, V_in, orbit_params);

    % Initial conditions
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0; y0; dxdt0; dydt0];

    h_ref = 0.000001;
    p = 4;
    % Time range and error tolerance values for the adaptive integrator
    t_range = [0, 30];
    error_range = logspace(-10, -1, 10);
    error_d = 0.001;
    % Set up Butcher tableau for Bogacki RK method
    Bogacki = struct();
    Bogacki.C = [0, 1/2, 3/4, 1];
    Bogacki.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
    Bogacki.A = [0, 0, 0, 0;
                 1/2, 0, 0, 0;
                 0, 3/4, 0, 0;
                 2/9, 1/3, 4/9, 0];

    % Preallocate arrays for storing results
    trunc_errors_adaptive = zeros(size(error_range));
    avg_step_sizes_adaptive = zeros(size(error_range));
    func_evals_adaptive = zeros(size(error_range));
    step_failure_rates = zeros(size(error_range));

    global_truncation_error_stepe( ...
            @(t, V) gravity_rate_func_tabby(t, V, orbit_params), ...
            @(t) compute_planetary_motion(t, V0, orbit_params), t_range, V0, Bogacki, 4, error_range, 'adaptive');
    % Adaptive step size integration over error tolerances
    for i = 1:length(error_range)
        % Run adaptive step integrator
        [t_list_adaptive, V_list_adaptive, avg_step_adaptive, num_evals_adaptive, fail_rate] = ...
            explicit_RK_variable_step_integration(gravity_rate_func, t_range, V0, h_ref, Bogacki, p, error_range(i));

        % Compute truncation error for adaptive step
        trunc_errors_adaptive(i) = global_truncation_error_stepe( ...
            @(t, V) gravity_rate_func_tabby(t, V, orbit_params), ...
            @(t) compute_planetary_motion(t, V0, orbit_params), t_range, V0, Bogacki, 4, error_range(i), 'adaptive');
        avg_step_sizes_adaptive(i) = avg_step_adaptive;
        func_evals_adaptive(i) = num_evals_adaptive;
        step_failure_rates(i) = fail_rate;
    end
    %Fixed step integrator over range of step sizes
    h_values = logspace(-5, -1, 10);
    trunc_errors_fixed = zeros(size(h_values));
    avg_step_sizes_fixed = h_values;  % For fixed step, avg step size is simply h
    func_evals_fixed = zeros(size(h_values));

    for j = 1:length(h_values)
        [t_list_fixed, V_list_fixed, ~, num_evals_fixed] = ...
            fixed_step_integration_stepe(gravity_rate_func, @(f, t, x, h) explicit_RK_step_stepe(f, t, x, h, Bogacki), t_range, V0, h_values(j));

        % Compute truncation error for fixed step
        trunc_errors_fixed(j) = global_truncation_error_stepe( ...
            @(t, V) gravity_rate_func_tabby(t, V, orbit_params), ...
            @(t) compute_planetary_motion(t, V0, orbit_params), t_range, V0, Bogacki, 4, error_d, 'fixed');
        func_evals_fixed(j) = num_evals_fixed;
    end

    % Plot truncation error vs. avg step size
    figure;
    loglog(avg_step_sizes_adaptive, trunc_errors_adaptive, 'bo-', 'LineWidth', 2, 'DisplayName', 'Adaptive Step');
    hold on;
    loglog(avg_step_sizes_fixed, trunc_errors_fixed, 'rx-', 'LineWidth', 2, 'DisplayName', 'Fixed Step');
    xlabel('Average Step Size (h)');
    ylabel('Global Truncation Error');
    legend('Location', 'Best');
    title('Global Truncation Error vs. Step Size for Adaptive and Fixed Step Integrators');
    grid on;
    hold off;

    % Plot global truncation error vs. number of function evaluations
    figure;
    loglog(func_evals_adaptive, trunc_errors_adaptive, 'bo-', 'LineWidth', 2, 'DisplayName', 'Adaptive Step');
    hold on;
    loglog(func_evals_fixed, trunc_errors_fixed, 'rx-', 'LineWidth', 2, 'DisplayName', 'Fixed Step');
    xlabel('Number of Function Evaluations');
    ylabel('Global Truncation Error');
    legend('Location', 'Best');
    title('Global Truncation Error vs. Function Evaluations for Adaptive and Fixed Step Integrators');
    grid on;
    hold off;
    
    % Plot failure rate vs. average step size for adaptive step integrator
    figure;
    semilogx(avg_step_sizes_adaptive, step_failure_rates, 'bo-', 'LineWidth', 2);
    xlabel('Average Step Size (h)');
    ylabel('Failure Rate');
    title('Failure Rate vs. Average Step Size for Adaptive Step Integrator');
    grid on;

end
