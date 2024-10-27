%Runs numerical integration arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration


function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration_tabby ...
(rate_func_in,tspan,X0,h_ref,BT_struct)

    % calculate steps and h
    [num_steps, h_avg] = iteration_solver(tspan, h_ref);
    
    % define variables
    XA = X0
    num_evals = 0;
    t_list = linspace(tspan(1),tspan(2),num_steps+1);
    X_list = zeros(num_steps+1,length(X0))
    X_list(1,:) = X0';
    % calculate the values until it is just short of the end value
    for i = 1:num_steps
        t = t_list(i);
        % [XB, temp_eval] = explicit_RK_step_tabby(rate_func_in,t,XA,h_avg, BT_struct);
        [XB1, XB2, temp_eval] = RK_step_embedded(rate_func_in,t,XA,h_avg,BT_struct)
        num_evals = num_evals + temp_eval;
        X_list(i,:) = XB1;
        X_list(i+1,:) = XB2;
    end
end


function [num_steps, h] = iteration_solver(tspan, h_ref)
    range = tspan(2)-tspan(1);
    num_steps = range/h_ref; % The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps); % Round the number of steps up (to get a real number)
    h = range/num_steps; % Divide range by steps to get real h
end