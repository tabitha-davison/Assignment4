function [XB, num_evals] = explicit_RK_step_stepe(rate_func_in,t,XA,h,BT_struct)

    s = length(BT_struct.C);               
    K = zeros(length(XA), s);            % Stage evaluations
    num_evals = 0;                       % Function evaluation counter
    
    % Loop through each stage to calculate K values
    for i = 1:s
        % Compute sum of previous stages for intermediate state calculation
        stage_sum = XA;
        for j = 1:i-1
            stage_sum = stage_sum + h * BT_struct.A(i, j) * K(:, j);
        end
        
        % Intermediate time and state for this stage
        t_i = t + BT_struct.C(i) * h;
        K(:, i) = rate_func_in(t_i, stage_sum);  % Compute derivative at stage
        num_evals = num_evals + 1;
    end
    
    % Calculate XB (next state) using the weighted sum of stages
    XB = XA;
    for i = 1:s
        XB = XB + h * BT_struct.B(i) * K(:, i);
    end
end
    

