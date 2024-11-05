function [t_list, X_list, h_avg, num_evals] = ...
    fixed_step_integration_stepe(rate_func_in, step_func, tspan, X0, h_ref)

    % Unpack tspan (start and end of the time interval)
    t_start = tspan(1);
    t_end = tspan(2);

    % Calculate the number of steps N and actual step size h_avg
    N = ceil((t_end - t_start) / h_ref);  % Number of integration steps
    h_avg = (t_end - t_start) / N;        % Actual time step used

    % Initialize arrays to store time values (t_list) and X values (X_list)
    t_list = linspace(t_start, t_end, N + 1)';  % Array of time steps
    X_list = zeros(N + 1, length(X0));          % Array to store X values

    % Set the initial condition (X at t_start)
    X_list(1, :) = X0';  % Ensure X0 is stored as a row vector

    % Initialize evaluation counter
    num_evals = 0;

    % Loop through each time step
    for i = 1:N
        % Compute the next step using the provided step function
        [X_next, evals] = step_func(rate_func_in, ...
            t_list(i), X_list(i, :)', h_avg);

        % Store the result for the next time step
        X_list(i + 1, :) = X_next';
        
        % Accumulate the total number of function evaluations
        num_evals = num_evals + evals;
    end

end

