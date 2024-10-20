

BT_struct = struct();
BT_struct.A = [0, 0; 0.5, 0]; % A matrix
BT_struct.B = [0; 1]; % B vector
BT_struct.C = [0; 0.5]; % C vector

BT_struct


orbit_params = struct();

orbit_params.m_sun = 1.9*10^30;
orbit_params.m_planet = 5.9*10^24;
orbit_params.G = 9.81;

orbit_params