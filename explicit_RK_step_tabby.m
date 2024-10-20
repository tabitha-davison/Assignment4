

%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step

function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)

    % Get the number of stages in the Runge-Kutta method (s is length of C or rows of A)
    s = length(BT_struct.C);
    
    % Initialize the K matrix (each column represents a stage's evaluation)
    K = zeros(length(XA), s);

    for i = 1:s
        sum_ki = K * BT_struct.A(i, :)'; % Weighted sum of previous ks
        K(:, i) = rate_func_in(t + BT_struct.C(i) * h, XA + h * sum_ki);
    end
    
    % Compute the next value using weighted sum of ks
    XB = XA + h * (K * BT_struct.B);
    
    num_evals = 1;

end


