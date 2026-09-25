%% Load parameters (these two functions are treated as fixed inputs)
env  = get_env_params();
moth = get_boat_params();

% moth.b_sail is NaN in get_boat_params.m -- set an assumed rig span
% here until you have (or are allowed to add) the real value:
moth.b_sail = 3.5;    % [m] <-- PLACEHOLDER, replace with real sail span

%% Define analysis condition
tw_deg    = 10;               % [deg] total twist
tw        = tw_deg*pi/180;    % [rad]
VS        = 15;               % [m/s] boat speed
tws       = 8;                % [m/s] true wind speed
twa_deg   = 60;                % [deg] true wind angle
twa       = twa_deg*pi/180;   % [rad]
phi_deg   = 0;                 % [deg] heel angle
phi       = phi_deg*pi/180;   % [rad]
theta_deg = 0;                 % [deg] pitch angle
theta     = theta_deg*pi/180; % [rad]

%% Run
[FA, CEA, alfa] = calc_sail(moth, env, tw, VS, tws, twa, phi, theta);

fprintf('F_D (driving)   = %.1f N\n', FA(1));
fprintf('F_SF (side)     = %.1f N\n', FA(2));
fprintf('F_vert          = %.1f N\n', FA(3));
fprintf('CE height       = %.2f m above deck\n', CEA);
fprintf('alpha           = %.2f deg\n', alfa*180/pi);
