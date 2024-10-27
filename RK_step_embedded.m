function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)

%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
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
%XB1: the approximate value for X(t+h) using the first row of the Tableau
%XB2: the approximate value for X(t+h) using the second row of the Tableau
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step

    % Initialize the K matrix (each column represents a stage's evaluation)
    K = zeros(length(XA), length(BT_struct.C));
    K(:,1) = rate_func_in(t,XA);
    
    for i = 2:length(BT_struct.C)
        % K
        % BT_struct.A(i, :)'
        sum_ki = K * BT_struct.A(i, :)'; % Weighted sum of previous ks
        K(:, i) = rate_func_in(t + (BT_struct.C(i) * h), XA + (h * sum_ki));
    end
    
    % Compute the next value using weighted sum of ks
    XB1 = XA + h * (K * BT_struct.B(1,:)');
    XB2 = XA + h * (K * BT_struct.B(2,:)');
    num_evals = length(BT_struct.B);

end


