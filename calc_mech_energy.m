function E = calc_mech_energy(V, params)
    x = V(1);
    y = V(2);
    vx = V(3);
    vy = V(4);
    
    KE = 0.5 * params.m_planet * (vx^2 + vy^2);
    
    r = sqrt(x^2 + y^2);
    PE = -params.G * params.m_sun * params.m_planet / r;
    
    % Total mechanical energy
    E = KE + PE;
end