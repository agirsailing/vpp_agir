function [FA, CEA, alfa_vec] = calc_sail(moth, env, tw, VS, tws, twa, phi, theta)
% CALC_SAIL  Aerodynamic sail force in the boat reference frame.
%
%   [FA, CEA, alfa_vec] = calc_sail(moth, env, tw, VS, tws, twa, phi, theta)
%
% INPUTS
%   moth  [struct] boat/rig geometry, from get_boat_params()
%   env   [struct] environment (air density, wind), from get_env_params()
%   tw    [rad]    total sail twist, root-to-tip
%   VS    [m/s]    boat speed
%   tws   [m/s]    true wind speed
%   twa   [rad]    true wind angle (0 = dead upwind)
%   phi   [rad]    heel angle, positive = heeling to leeward   (optional, default 0)
%   theta [rad]    pitch angle, positive = bow up              (optional, default 0)
%
% OUTPUTS
%   FA       [N]   3x1 force in boat axes [F_D; F_SF; F_vert]
%                  (x fwd, y stbd, z down -- matches get_boat_params.m convention)
%   CEA      [m]   height of the sail centre of effort above the deck (+up)
%   alfa_vec [rad] angle of attack used in the CL/CD model (currently a
%                  single mean value -- kept as a "vec" so this can be
%                  extended to a spanwise strip-theory loop later)
%
% NOTES / THINGS YOU STILL NEED TO CONFIRM OR TUNE
%   - moth.length_sail is NaN in get_boat_params.m (which you said can't be
%     edited), so the caller MUST set an assumed sail span, e.g.
%         moth.length_sail = 3.5;   % [m] <-- replace with the real rig span
%     before calling this function, otherwise it errors out on purpose
%     rather than silently returning NaN.
%   - e0 (span efficiency), flat (flattening) and CD0 (parasitic drag)
%     below are placeholders -- tune them to your sail's real polars.
%   - The angle-of-attack model (alpha = awa_eff - tw/2) is a crude
%     "average station" stand-in for a proper strip-theory integration
%     over the twisted span. Good enough for a first pass, not for a
%     final VPP.
%   - asin() in the velocity triangle only returns the right AWA up to
%     90 deg; see the commented atan2 alternative below if you need to
%     run downwind cases (broad TWA).

if nargin < 7 || isempty(phi),   phi   = 0; end
if nargin < 8 || isempty(theta), theta = 0; end

%% 1) Apparent wind from the velocity triangle (horizontal boat plane)
aws = sqrt(tws^2 + VS^2 - 2*tws*VS*cos(pi - twa));   % [m/s]
awa = asin( sin(pi - twa) * tws / aws );              % [rad] -- valid for AWA <= 90 deg

% For AWA that can exceed 90 deg (reaching/running), use instead:
% Vx  = tws*cos(twa) + VS;
% Vy  = tws*sin(twa);
% aws = hypot(Vx, Vy);
% awa = atan2(Vy, Vx);

%% 2) Effective wind seen by the (heeled, pitched) rig -- eq. (34)-(35)
awa_eff = atan( (tan(awa)*cos(phi) - sin(theta)*sin(phi)) / cos(theta) );
V_eff   = aws * sqrt( (sin(awa)*cos(phi) - cos(awa)*sin(phi)*sin(theta))^2 ...
                     + cos(awa)^2 * cos(theta)^2 );

%% 3) Angle of attack and aero coefficients (evaluated with the EFFECTIVE wind)
alpha    = awa_eff - tw/2;             % crude mean-station AoA, see notes above
alfa_vec = alpha;

rho_air = env.air.rho;                 % [kg/m^3]
q_eff   = 0.5 * rho_air * V_eff^2;     % [Pa]

A = moth.S_sail;                       % [m^2]
if isnan(moth.length_sail)
    error(['moth.length_sail (sail span) is NaN -- get_boat_params.m does not ' ...
           'set it. Assign an assumed span to moth.length_sail before calling ' ...
           'calc_sail, e.g. moth.length_sail = 3.5; %% [m]']);
end
AR = moth.length_sail^2 / A;                % [-]

a0   = 2*pi;                           % [1/rad] thin-aerofoil 2D lift slope
e0   = 0.85;                           % [-] span efficiency  -- TUNE
flat = 1.0;                            % [-] flattening (1 = full power) -- TUNE
CLa  = a0 / (1 + a0/(pi*e0*AR));       % [1/rad] 3D lift slope
CL   = flat * CLa * alpha;             % [-]
if ~isnan(moth.CLmax_sail)
    CL = min(CL, moth.CLmax_sail);     % crude stall cap
end
CD0  = 0.008;                          % [-] parasitic drag -- TUNE
CDi  = CL^2 / (pi*e0*AR);              % [-] induced drag
CD   = CD0 + CDi;                      % [-]

L = q_eff * A * CL;                    % [N] lift, perpendicular to apparent wind
D = q_eff * A * CD;                    % [N] drag, along apparent wind

%% 4) Resolve into boat axes
% Standard VPP decomposition (matches F_D / F_SF in your diagram):
% rotate lift/drag by the ACTUAL (not effective) AWA to land in the
% horizontal boat x-y plane.
F_D  = L*sin(awa) - D*cos(awa);        % [N] driving force, +fwd
F_SF = L*cos(awa) + D*sin(awa);        % [N] side (heeling) force, + to leeward

% Heel and pitch tilt that horizontal-plane force out of the x-y plane.
R_heel  = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
R_pitch = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
FA = R_pitch * R_heel * [F_D; F_SF; 0];

%% 5) Centre of effort
CEA = -moth.z_CE_sail;                 % [m] height above deck (z is +down in moth frame)

end
