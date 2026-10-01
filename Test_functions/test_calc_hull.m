% Test of calc_hull with a table from make_preCalcHull
clear; close all
env = get_env_params(); boat = get_boat_params();
rho = env.water.rho; g = env.g; m = boat.m_tot;
hulldata.preCalc = 'preCalcHull_kayak.mat';
load(hulldata.preCalc)

% Static float (should be minus the draft from make_preCalcHull)
Hs = fzero(@(H) rho*hull_V(H, hulldata) - m, [H_(1) -1e-6]);
fprintf('Static float: H = %.4f m\n', Hs);

% Limits (should all be 0)
[R0, V0] = calc_hull(3, 0, hulldata);
[R1, V1] = calc_hull(3, 0.1, hulldata);
fprintf('H = 0: R = %g, V = %g | H = 0.1: R = %g, V = %g\n', R0, V0, R1, V1);

% Take-off sweep: foils lift W*(VS/VS_to)^2, hull carries the rest
VS_to = 4;
VS = linspace(0.5, 5, 60);
H = zeros(size(VS)); R = zeros(size(VS));
for j = 1:numel(VS)
    Wh = m*g*max(1 - (VS(j)/VS_to)^2, 0);
    if Wh > 0
        H(j) = fzero(@(h) rho*g*hull_V(h, hulldata) - Wh, [H_(1) -1e-9]);
        R(j) = calc_hull(VS(j), H(j), hulldata);
    end
end
subplot(2,1,1); plot(VS, H, LineWidth=1.2); grid on; ylabel('H [m]')
subplot(2,1,2); plot(VS, R, LineWidth=1.2); grid on; ylabel('R hull [N]'); xlabel('VS [m/s]')

function V = hull_V(H, hulldata)
    [~, V] = calc_hull(1, H, hulldata);
end