function make_preCalcHull(brifile, outfile, plotflag, z_down, bow_at_max_x)
% Hull table (R_, V_, H_, VS_) for calc_hull.m, Findlay & Turnock (2008):
% ITTC-57 friction + Delft (Keuning & Katgert 2008) residuary resistance.
% Hydrostatics with the course functions (ReadHullGeometry, WetSections,
% CalculateHydrostatics). H < 0 = hull bottom below the water.

if nargin < 3, plotflag = false; end
if nargin < 4, z_down = true; end
if nargin < 5, bow_at_max_x = true; end

thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot);
addpath(genpath(thisDir));


env = get_env_params();
boat = get_boat_params();
rho = env.water.rho;
nu = env.water.nu;
g = env.g;
m = boat.m_tot;
k_hull = 0;
VS_ = 0.5:0.25:12;

%% Geometry
ship = ReadHullGeometry(brifile);
if z_down
    ship.z = -ship.z; % course tools use z up
end

ship.KG = boat.z_boat_cg;
ship.zk = zeros(1, ship.ns);
zmax = -inf;

for k = 1:ship.ns
    ship.zk(k) = min(ship.z(1:ship.np(k), k));
    zmax = max(zmax, max(ship.z(1:ship.np(k), k)));
end

ship.zkeel = min(ship.zk);
depth = zmax - ship.zkeel;
x = ship.x(1,:);
fprintf('Hull read: %d sections, x = %.3f to %.3f m, max beam %.3f m, depth %.3f m\n', ...
    ship.ns, min(x), max(x), 2*max(abs(ship.y(:))), depth);

if max(x) - min(x) > 10
    warning('Hull longer than 10 m: is the .bri in mm?');
end

if z_down && abs(ship.zkeel + boat.z_hull_bottom) > 0.05
    warning('Keel at z = %.3f m, expected %.3f m. Is the .bri really z down?', ...
        -ship.zkeel, boat.z_hull_bottom);
end

%% Draft at full weight
Vtarget = m/rho;
dscan = linspace(0, 0.999*depth, 41); % WetSections fails if the WL hits offset points exactly
Vscan = zeros(size(dscan));

for i = 2:numel(dscan)
    Vscan(i) = hydro(ship, -dscan(i), env, bow_at_max_x).V;
end

i2 = find(Vscan >= Vtarget, 1);

if isempty(i2)
    fprintf('  %8.4f   %9.5f\n', [dscan; Vscan]);
    error('Hull cannot carry %.1f kg: max volume %.4f m3, needed %.4f m3.', m, max(Vscan), Vtarget);
end

d0 = fzero(@(d) hydro(ship, -d, env, bow_at_max_x).V - Vtarget, dscan([i2-1 i2]));
full = hydro(ship, -d0, env, bow_at_max_x);

%% Table
H_ = linspace(-1.25*d0, 0, 26);
nH = numel(H_);
V_ = zeros(1, nH); S_ = zeros(1, nH); L_ = zeros(1, nH);

for i = 1:nH
    hs = hydro(ship, H_(i), env, bow_at_max_x);
    V_(i) = hs.V; S_(i) = hs.S; L_(i) = hs.Lwl;
end

[RrPerW, Fn_tab] = delft_rr_per_weight(full);
Fn = VS_ / sqrt(g*full.Lwl);
Fn_hold = 0.70; % Fn = 0.75 row gives Rr ~ 0, not used
keep = Fn_tab <= Fn_hold + 1e-9;
rr = interp1([0 Fn_tab(keep)], [0 RrPerW(keep)], min(Fn, Fn_hold), 'pchip');

R_ = zeros(nH, numel(VS_));

for i = 1:nH
    if V_(i) <= 0, continue; end
    Re = VS_ * 0.7*L_(i) / nu;
    Cf = 0.075 ./ (log10(Re) - 2).^2;
    Rf = 0.5*rho*VS_.^2 * S_(i).*Cf*(1+k_hull);
    Rr = rr * rho*g*V_(i);
    R_(i,:) = Rf + Rr;
end

save(outfile, 'R_', 'V_', 'H_', 'VS_');

fprintf('Mass %.1f kg -> draft %.3f m, LWL %.3f m, BWL %.3f m\n', m, d0, full.Lwl, full.Bwl);
fprintf('Swet %.3f m2, Awp %.3f m2, Cp %.3f, Cm %.3f, LCB/LWL %.3f\n', ...
    full.S, full.Aw, full.Cp, full.Cm, full.LCBfpp/full.Lwl);

if plotflag
    figure; subplot(1,2,1); plot(-H_, V_*rho, 'o-', LineWidth=1.2); grid on
    xlabel('Immersion -H [m]'); ylabel('Displacement [kg]');
    subplot(1,2,2); plot(VS_, R_(1:5:end,:)', LineWidth=1.2); grid on
    xlabel('VS [m/s]'); ylabel('R [N]');
    legend(arrayfun(@(h) sprintf('H=%.3f', h), H_(1:5:end), 'UniformOutput', false));
    state.eta = [0 0 -d0-ship.zkeel 0 0 0]; state.plotflag = 1;
    WetSections(ship, state);
end
end

function hs = hydro(ship, H, env, bow_at_max_x)
state.eta = [0 0 H-ship.zkeel 0 0 0]; % do not set state.H (wave height)
state.rau = env.water.rho;
state.g = env.g;
HS = CalculateHydrostatics(ship, state);
WS = WetSections(ship, state);
x = ship.x(1,:);

hs.V = abs(HS.V);
hs.S = trapz(x, WS.Ss);
hs.Aw = trapz(x, WS.bwl);
hs.Bwl = max(WS.bwl);
hs.Ax = max(abs(WS.Aws));
hs.Tc = max(-H, 0);
hs.Cm = hs.Ax / max(hs.Bwl*hs.Tc, eps);
hs.Lwl = 0; hs.Cp = NaN; hs.LCBfpp = NaN; hs.LCFfpp = NaN;

zwl = ship.zkeel - H;
zk = ship.zk; n = ship.ns;
wet = zk < zwl;
if hs.V <= 0 || ~any(wet)
    hs.V = 0; hs.S = 0; hs.Aw = 0; hs.Bwl = 0; hs.Ax = 0;
    return
end
i1 = find(wet, 1, 'first'); i2 = find(wet, 1, 'last');

xa = x(i1); xf = x(i2);

if i1 > 1, xa = interp1(zk([i1-1 i1]), x([i1-1 i1]), zwl); end
if i2 < n, xf = interp1(zk([i2 i2+1]), x([i2 i2+1]), zwl); end
hs.Lwl = xf - xa;
hs.Cp = hs.V / (hs.Lwl*hs.Ax);
xB = HS.CoB(1);
xF = trapz(x, x.*WS.bwl) / hs.Aw;

if bow_at_max_x
    hs.LCBfpp = xf - xB; hs.LCFfpp = xf - xF;
else
    hs.LCBfpp = xB - xa; hs.LCFfpp = xF - xa;
end
end

function [RrPerW, Fn] = delft_rr_per_weight(h)
% Rr/(rho g V) = (a0 + a1*LCB/L + a2*Cp + a3*V^(2/3)/Aw + a4*B/L
%                + a5*LCB/LCF + a6*B/Tc + a7*Cm) * V^(1/3)/L
[Fn, a] = delft_coefficients();
Fn = Fn.';
L = h.Lwl; V = h.V; s = V^(1/3)/L;
X = [1, h.LCBfpp/L, h.Cp, V^(2/3)/h.Aw, h.Bwl/L, h.LCBfpp/h.LCFfpp, h.Bwl/h.Tc, h.Cm];

names = {'LCB/LWL','Cp','V^(2/3)/Aw','BWL/LWL','LCB/LCF','BWL/Tc','Cm','V^(1/3)/LWL'};
val = [X(2:8) s];
lo = [0.500 0.521 0.079 0.170 0.930 2.46 0.646 0.120];
hi = [0.579 0.580 0.265 0.366 1.002 19.38 0.790 0.230];
fprintf('Delft applicability (value / series range):\n');
for i = 1:numel(val)
    flag = ''; if val(i) < lo(i) || val(i) > hi(i), flag = '  <-- OUTSIDE'; end
    fprintf('  %-12s %7.3f   [%6.3f %6.3f]%s\n', names{i}, val(i), lo(i), hi(i), flag);
end

RrPerW = max(((a*X.') * s).', 0);
end

function [Fn, a] = delft_coefficients()
% Keuning & Katgert (2008), Table 2. Rows: Fn; columns: a0...a7.
Fn = (15:5:75)' / 100;
a = [ ...
    -0.0005  0.0023 -0.0086 -0.0015  0.0061  0.0010  0.0001  0.0052;
    -0.0003  0.0059 -0.0064  0.0070  0.0014  0.0013  0.0005 -0.0020;
    -0.0002 -0.0156  0.0031 -0.0021 -0.0070  0.0148  0.0010 -0.0043;
    -0.0009  0.0016  0.0337 -0.0285 -0.0367  0.0218  0.0015 -0.0172;
    -0.0026 -0.0567  0.0446 -0.1091 -0.0707  0.0914  0.0021 -0.0078;
    -0.0064 -0.4034 -0.1250  0.0273 -0.1341  0.3578  0.0045  0.1115;
    -0.0218 -0.5261 -0.2945  0.2485 -0.2428  0.6293  0.0081  0.2086;
    -0.0388 -0.5986 -0.3038  0.6033 -0.0430  0.8332  0.0106  0.1336;
    -0.0347 -0.4764 -0.2361  0.8726  0.4219  0.8990  0.0096 -0.2272;
    -0.0361  0.0037 -0.2960  0.9661  0.6123  0.7534  0.0100 -0.3352;
     0.0008  0.3728 -0.3667  1.3957  1.0343  0.3230  0.0072 -0.4632;
     0.0108 -0.1238 -0.2026  1.1282  1.1836  0.4973  0.0038 -0.4477;
     0.1023  0.7726  0.5040  1.7867  2.1934 -1.5479 -0.0115 -0.0977];
end