function integration_assign4 ()
    BT_struct.A = [0 0; 0.5 0];
    BT_struct.B = [0; 1];
    BT_struct.C = [0; 0.5];
    tspan = [0, 20];
    X0 = 1;
    t = linspace(0, 20, 100000);
    h_ref = 0.38;

[t_list, X_list, h_avg, num_evals] = ...
    fixed_step_integration(@rate_func01, @(f,t,X,h) explicit_RK_step(f,t,X,h,BT_struct), tspan, X0, h_ref)
h_avg
end

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end
% % function to compute the value of X at next time step
% % using explicit midpoint approximation
% % 
% % INPUTS:
% % rate_func_in: the function used to compute dXdt. rate_func_in will
% %  have the form: dXdt = rate_func_in(t,X) (t is before X)
% % t: the value of time at the current step
% % XA: the value of X(t)
% % h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
% % 
% % OUTPUTS:
% % XB: the approximate value for X(t+h) (the next step)
% %  formula depends on the integration method used
% % num_evals: A count of the number of times that you called
% %  rate_func_in when computing the next step
% 
% function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
% 
%     % compute midpoint approximation
%     % Xn+0.5 = Xn + (h/2) * f(tn, Xn)
%     half_step = h / 2;
%     dXdt_n = rate_func_in(t, XA);
%     X_mid = XA + half_step * dXdt_n;
% 
%     % compute Xn +  using midpoint approximation
%     % Xn+1 = Xn + h * f(tn + h/2, Xn+0.5)
%     dXdt_mid = rate_func_in(t + half_step, X_mid);
%     XB = XA + h * dXdt_mid;
% 
%     % add 2 to num_evals
%     num_evals = 2;
% 
% end

% function to compute the value of X at the next time step
% using any explicit Runge-Kutta method (including midpoint, Euler, etc.)
% specified by the Butcher tableau.
% 
% INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
%   have the form: dXdt = rate_func_in(t,X) (t is before X)
% t: the value of time at the current step
% XA: the value of X(t)
% h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
% BT_struct: a struct that contains the Butcher tableau
%   BT_struct.A: matrix of a_{ij} values
%   BT_struct.B: vector of b_i values
%   BT_struct.C: vector of c_i values
% OUTPUTS:
% XB: the approximate value for X(t+h) (the next step)
% num_evals: A count of the number of times that you called
%   rate_func_in when computing the next step

function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)

    % Get the number of stages in the Runge-Kutta method (s is length of C or rows of A)
    s = length(BT_struct.C);
    
    % Initialize the K matrix (each column represents a stage's evaluation)
    K = zeros(length(XA), s);
    
    % Number of function evaluations (to keep track)
    num_evals = 0;
    
    % Loop over the number of stages
    for i = 1:s
        % Compute the sum of a_{ij} * k_j from the previous stages
        sum_val = 0;
        for j = 1:i-1
            sum_val = sum_val + BT_struct.A(i,j) * K(:,j);
        end
        
        % Compute the argument for f(t + c_i*h, X + sum_val)
        t_i = t + BT_struct.C(i) * h;
        X_i = XA + h * sum_val;
        
        % Evaluate the derivative f(t_i, X_i)
        K(:,i) = rate_func_in(t_i, X_i);
        
        % Increment the function evaluation count
        num_evals = num_evals + 1;
    end
    
    % Compute the next step using the weights b_i
    XB = XA + h * K * BT_struct.B(:);

end


function [t_list, X_list, h_avg, num_evals] = ...
    fixed_step_integration(rate_func_in, step_func, tspan, X0, h_ref)

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

