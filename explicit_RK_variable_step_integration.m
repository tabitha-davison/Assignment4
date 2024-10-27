
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration

function [t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct, p, error_desired)

    % calculate steps and h
    [num_steps, h_next] = iteration_solver(tspan, h_ref);
    
    % define variables
    XA = X0
    num_evals = 0;
    t_current = 0;
    t_list = tspan(1);
    X_list = X0;
    % t_list = linspace(tspan(1),tspan(2),num_steps+1);
    % X_list = zeros(num_steps+1,length(X0))
    % X_list(1,:) = X0';

    % calculate the values until it is just short of the end value
    while t_current < tspan(2)
        t_next = t_current + h_next;

        if t_next > tspan(2)
            h_next = tspan(2)- t_current;
            t_next = tspan(2);
        end

        % [XB, temp_eval] = explicit_RK_step_tabby(rate_func_in,t,XA,h_avg, BT_struct);
        [XB, evals, h_next, redo] = explicit_RK_variable_step...
         (rate_func_in, t_next, XA, h_next, BT_struct, p, error_desired)
        if ~ redo
            t_current = t_next;
            XA = XB;
            X_list(:, end+1) = XB;
            t_list(end+1) = t_current;
            num_evals = num_evals + evals;
        end
    end
    h_avg = (t_list(end) - t_list(1))/(length(t_list) - 1);
end


function [num_steps, h] = iteration_solver(tspan, h_ref)
    range = tspan(2)-tspan(1);
    num_steps = range/h_ref; % The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps); % Round the number of steps up (to get a real number)
    h = range/num_steps; % Divide range by steps to get real h
end