function [R, V, xB] = calc_hull(VS, H, hulldata)
% Resistance, displaced volume and centre of buoyancy of the hull,
% interpolated from the table made by make_preCalcHull.
% H = 0 hull out of the water, H < 0 immersed.
% H below the deepest table row is clamped to that row (never returns 0).
%
% The table is loaded once and kept in memory (persistent). If you
% regenerate the .mat file under the same name, run "clear calc_hull".
%
% Outputs
% R   [N]   resistance
% V   [m3]  displaced volume
% xB  [m]   x of the centre of buoyancy (boat coordinates)

persistent tab file

if isempty(tab) || ~strcmp(file, hulldata.preCalc)      % First call or new file
    d = load(hulldata.preCalc);                         % R_ (rows H_, columns VS_), V_, xB_, H_, VS_
    tab.R    = griddedInterpolant({d.H_, d.VS_}, d.R_, 'linear', 'linear');
    tab.V    = griddedInterpolant(d.H_, d.V_, 'linear');
    tab.xB   = griddedInterpolant(d.H_, d.xB_, 'linear');
    tab.Hmin = d.H_(1);
    file = hulldata.preCalc;
end

if H >= 0                       % Hull out of the water
    R = 0; V = 0; xB = 0;
    return
end

if H < tab.Hmin                 % Deeper than the table: clamp, do not drop to zero
    warning('calc_hull:clamped', ...
        'H = %.3f m is below the table limit %.3f m, clamped.', H, tab.Hmin);
    H = tab.Hmin;
end

R  = max(tab.R(H, VS), 0);      % Linear extrapolation in VS, can go negative at low VS
V  = max(tab.V(H), 0);
xB = tab.xB(H);
end