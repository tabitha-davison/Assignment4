
BT_struct = struct();
BT_struct.A = [0, 0; 0.5, 0]; % A matrix
BT_struct.B = [0; 1]; % B vector
BT_struct.C = [0; 0.5]; % C vector

BT_struct

tspan = [0, 3.154e7];
h_ref = 86400; % 1 day (in seconds)
t_list = linspace(tspan(1),tspan(2),100+1);
XA = 1;
h = 0.38;

% function dVdt = gravity_rate_func(t,V,orbit_params)
% 
%     xp = V(1); % Position in x
%     yp = V(2); % Position in y
%     vx = V(3); % Velocity in x
%     vy = V(4); % Velocity in y
% 
%     G = orbit_params.G; % Gravitational constant
%     m_sun = orbit_params.m_sun; % Mass of the sun (kg)
%     m_planet = orbit_params.m_planet; % Mass of the planet (kg)
% 
%     r = sqrt(xp^2 + yp^2);
%     
%     accel_factor = -G * m_sun / r^3; 
% 
%     ax = accel_factor * xp; % Acceleration in x direction
%     ay = accel_factor * yp; % Acceleration in y direction
% 
%     % Construct the derivative of the state vector dVdt
%     dVdt = [vx; vy; ax; ay];
% 
% end


[XB, num_evals] = explicit_RK_step(@gravity_rate_func,t_list,XA,h,BT_struct)















% 
% [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
% (gravity_rate_func_tabby,tspan,X0,h_ref,orbit_params)

