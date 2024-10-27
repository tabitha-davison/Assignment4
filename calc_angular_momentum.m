function H = calc_angular_momentum(V, params)
    x = V(1);
    y = V(2);
    vx = V(3);
    vy = V(4);
    
    H = params.m_planet * (x * vy - y * vx);
end