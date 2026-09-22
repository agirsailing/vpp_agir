function delft_resistance(env, boat, moth, V)
% delft_resistance - Calculate the resistance of a hull using Delft method

%   Formula used:
% R_hull(V) = R_F(V) + R_R,Delft(V)_; where V is velocity


%   Inputs:

%  environmentals
rho = env.rho; % [kg/m^3] water density
nu = env.water.nu;   % [m^2/s] kinematic viscosity
g = env.water.g;     % [m/s^2] gravity acceleration

% boat parameters
L = boat.hull.LWL; % [m] waterline length
S = boat.hull.Swet; % [m^2] wetted surface area


% Froude number
Fn = V / sqrt(g*L);

% Friction
Re = V * L / nu; % Reynolds number
Cf = 0.075 / (log10(Re) - 2)^2;

Rf = 0.5 * rho * V^2 * S * Cf;

% Delft residuary resistance
%Rr = delftResistance(boat, Fn);

Rr = 0;

R = Rf + Rr;



end
