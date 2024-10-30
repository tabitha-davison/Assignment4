function experimental_plot()
    % Set up parameters for sun and planet
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    
    % Define the gravity rate function
    gravity_rate_func = @(t_in, V_in) gravity_rate_func_tabby(t_in, V_in, orbit_params);
    
    % Initial conditions
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0; y0; dxdt0; dydt0];
    
    % Set up time range and step sizes
    t_range = [0, 30];
    h_values = logspace(-5, -1, 50);
    
    % Define the Butcher tableau for the Bogacki-Shampine method
    Bogacki = struct();
    Bogacki.C = [0, 1/2, 3/4, 1];
    Bogacki.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
    Bogacki.A = [0, 0, 0, 0;
                 1/2, 0, 0, 0;
                 0, 3/4, 0, 0;
                 2/9, 1/3, 4/9, 0];
    % Initialize arrays to store deviations for each h value
    energy_dev = zeros(length(h_values), 1);
    momentum_dev = zeros(length(h_values), 1);

    % Loop over each h value
    for i = 1:length(h_values)
        h = h_values(i);

        % Integrate with current step size
        [t_list, V_list] = fixed_step_integration_stepe(gravity_rate_func, ...
                                                        @(f, t, X, h) explicit_RK_step_stepe(f, t, X, h, Bogacki), ...
                                                        t_range, V0, h);

        % Calculate initial energy and momentum
        initial_energy = calc_mech_energy(V0, orbit_params);
        initial_momentum = calc_angular_momentum(V0, orbit_params);

        % Initialize arrays to store energy and momentum values
        E_vals = zeros(length(t_list), 1);
        M_vals = zeros(length(t_list), 1);

        % Calculate energy and momentum for each step
        for j = 1:length(t_list)
            E_vals(j) = calc_mech_energy(V_list(j, :)', orbit_params);
            M_vals(j) = calc_angular_momentum(V_list(j, :)', orbit_params);
        end

        % Compute final deviation from initial values
        energy_dev(i) = abs(E_vals(end) - initial_energy);
        momentum_dev(i) = abs(M_vals(end) - initial_momentum);
    end

    % Plotting
    figure;

    % Plot energy deviation
    subplot(2, 1, 1);
    loglog(h_values, energy_dev, 'b-o', 'LineWidth', 1.5);
    title('Mechanical Energy Error vs Step Size h');
    xlabel('Step Size (h)');
    ylabel('Mechanical Energy Error');
    grid on;

    % Plot momentum deviation
    subplot(2, 1, 2);
    loglog(h_values, momentum_dev, 'r-o', 'LineWidth', 1.5);
    title('Angular Momentum Error vs Step Size h');
    xlabel('Step Size (h)');
    ylabel('Angular Momentum Error');
    grid on;
end

