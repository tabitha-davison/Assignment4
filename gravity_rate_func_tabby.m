%Rate function describing Newton’s law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
3
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)

    xp = V(1); % Position in x
    yp = V(2); % Position in y
    vx = V(3); % Velocity in x
    vy = V(4); % Velocity in y

    G = orbit_params.G; % Gravitational constant
    m_sun = orbit_params.m_sun; % Mass of the sun (kg)
    m_planet = orbit_params.m_planet; % Mass of the planet (kg)

    r = sqrt(xp^2 + yp^2);
    
    accel_factor = -G * m_sun / r^3; 

    ax = accel_factor * xp; % Acceleration in x direction
    ay = accel_factor * yp; % Acceleration in y direction

    % Construct the derivative of the state vector dVdt
    dVdt = [vx; vy; ax; ay];


%     tspan = [0, 3.154e7];
%     h_ref = 86400; % 1 day (in seconds)
%     
%     [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration( ...
%         @gravity_rate_func, tspan, X0, h_ref, BT_struct, orbit_params);
%     
%     figure;
%     plot(X_list(:, 1), X_list(:, 2));
%     xlabel('x position (m)');
%     ylabel('y position (m)');
%     title('Planetary Orbit');
%     grid on;
%     axis equal;

end