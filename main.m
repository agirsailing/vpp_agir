


% Get preset parameters
env = get_env_params(); % environmental parameters
moth = get_boat_params(); % boat parameters

% Boat parameters
boat.hull.LWL   = moth.L_wl;      % waterline length [m]
boat.hull.Swet  = 1;      % wetted surface area [m^2]
boat.hull.disp  = 30;      % displacement [kg] - intial estimate



