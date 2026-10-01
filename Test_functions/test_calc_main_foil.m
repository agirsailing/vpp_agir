% Test of calc_mainfoil: lift and drag vs boat speed at fixed foil angles
clear; clc; close all

moth = get_boat_params();
env = get_env_params();

% Operating point
h = moth.depth_foil_flying;     % [m] ride height
theta = 0;                      % [rad] trim
FY = 0;                         % [N] no side force, foil only
VS_vec = 0.1:0.1:10;            % [m/s] boat speeds
aoa = [0 1 2 3 4 5 6];          % [deg] foil angles to compare

% Main foil share of the weight
W = moth.m_total*env.g;
L_target = W*moth.x2/(moth.x1 + moth.x2);

% Sweep speed for each foil angle
L = zeros(length(aoa), length(VS_vec));
D = zeros(length(aoa), length(VS_vec));

for j = 1:length(aoa)
    for i = 1:length(VS_vec)
        [F, ~, ~] = calc_mainfoil(VS_vec(i), h, theta, aoa(j)*pi/180, FY, moth, env);
        L(j,i) = -F(3);         % [N] z down: lift = -Fz
        D(j,i) = -F(1);         % [N] drag = -Fx
    end
end

% Plots
VS_kn = VS_vec*env.ms2kn;
names = arrayfun(@(a) sprintf('AoA = %g deg', a), aoa, 'UniformOutput', false);

figure
subplot(3,1,1)
plot(VS_kn, L, 'LineWidth', 1.5); hold on
yline(L_target, 'k--', 'Main foil share of weight');
ylabel('Lift [N]'); grid on; legend(names, 'Location', 'northwest')
title(sprintf('Main foil, h = %.2f m', h))

subplot(3,1,2)
plot(VS_kn, D, 'LineWidth', 1.5)
ylabel('Drag [N]'); grid on

subplot(3,1,3)
plot(VS_kn, L./D, 'LineWidth', 1.5)
xlabel('VS [kn]'); ylabel('L/D [-]'); grid on
