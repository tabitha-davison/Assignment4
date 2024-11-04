function [global_error, h_list] = global_truncation_error_stepe(rate_func, solution_func, tspan, V0, Butcher_tableau, p, range, mode)
    % Global truncation error calculator for both fixed and adaptive methods.
    % rate_func        - function handle for the rate of change (differential equation)
    % solution_func    - function handle for the exact/analytical solution
    % tspan            - [t0, tf] time range for integration
    % V0               - Initial condition vector
    % Butcher_tableau  - Butcher tableau for Runge-Kutta method
    % p                - Order of the method
    % range            - Array of step sizes (for fixed) or error tolerances (for adaptive)
    % mode             - 'adaptive' for adaptive integrator, 'fixed' for fixed step

    % Preallocate arrays for step sizes and global errors
    h_list = zeros(length(range), 1);
    global_error = zeros(length(range), 1);

    % Loop through the desired range based on mode
    for n = 1:length(range)
        if strcmp(mode, 'adaptive')
            % Adaptive Step: Current desired error tolerance
            desired_error = range(n);
            h_ref = 0.0001;  % Reference initial step size for adaptive integrator

            % Run adaptive integrator
            [t_list, V_list, h_avg, ~] = explicit_RK_variable_step_integration( ...
                rate_func, tspan, V0, h_ref, Butcher_tableau, p, desired_error);

            % Record average step size
            h_list(n) = h_avg;

        elseif strcmp(mode, 'fixed')
            % Fixed Step: Current step size
            h_fixed = range(n);

            % Run fixed-step integrator
            [t_list, V_list, ~, ~] = fixed_step_integration_stepe( ...
                rate_func, @(f, t, x, h) explicit_RK_step_stepe(f, t, x, h, Butcher_tableau), ...
                tspan, V0, h_fixed);

            % Record fixed step size
            h_list(n) = h_fixed;
        end

        % Calculate the exact solution at each time step
        V_solution = solution_func(t_list);

        % Calculate the global error at the final time step
        global_error(n) = norm(V_list(end, :)' - V_solution(end, :)');
    end
end
