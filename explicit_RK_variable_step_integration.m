function [t_list, X_list, h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration ...
(rate_func_in, tspan, X0, h_ref, BT_struct, p, error_desired)

    % calculate steps and h
    [num_steps, h_next] = iteration_solver(tspan, h_ref);
    
    % define variables
    XA = X0;
    num_evals = 0;
    t_current = tspan(1);
    t_list = tspan(1);
    X_list = X0';
    
    % Counters for attempted and failed steps
    num_attempted_steps = 0;
    num_failed_steps = 0;

    % calculate the values until it is just short of the end value
    while t_current < tspan(2)
        t_next = t_current + h_next;

        if t_next > tspan(2)
            h_next = tspan(2) - t_current;
            t_next = tspan(2);
        end

        [XB, evals, h_next, redo] = explicit_RK_variable_step...
         (rate_func_in, t_next, XA, h_next, BT_struct, p, error_desired);

        % Increment attempted steps
        num_attempted_steps = num_attempted_steps + 1;
        
        if redo
            % Increment failed steps if redo is required
            num_failed_steps = num_failed_steps + 1;
        else
            % Successful step: update variables
            t_current = t_next;
            XA = XB;
            X_list(end+1,:) = XB';
            t_list(end+1) = t_current;
            num_evals = num_evals + evals;
        end
    end

    % Calculate average step size
    h_avg = (t_list(end) - t_list(1)) / (length(t_list) - 1);

    % Calculate step failure rate
    step_failure_rate = num_failed_steps / num_attempted_steps;
end


function [num_steps, h] = iteration_solver(tspan, h_ref)
    range = tspan(2) - tspan(1);
    num_steps = ceil(range / h_ref); % Round up to get a real number of steps
    h = range / num_steps; % Calculate actual h
end
