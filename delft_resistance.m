function [R, details] = delft_resistance(V, hull, env)
%DELFT_RESISTANCE Upright bare-hull resistance, Keuning-Katgert (2008).
%   [R, DETAILS] = DELFT_RESISTANCE(V, HULL, ENV) returns resistance in N
%   for a scalar or vector of speeds through water V [m/s]. Output arrays
%   retain V's shape. Nonzero speeds require 0.15 <= Fn <= 0.75.
%
%   HULL fields (scalar, SI): LWL, BWL, draft, volume, Swet, Awp, Cp, Cm,
%   xLCB, xLCF, xFPP. Coordinates are x-forward from the stern on deck;
%   xFPP is the forward static-waterline endpoint in that same frame.
%   ENV fields: water.rho [kg/m^3], water.nu [m^2/s], g [m/s^2].
%
%   DETAILS contains Rf, Rr, Fn, Re, Cf, model and validity. At V=0 the
%   forces and Re are zero and Cf is NaN (not evaluated). Rr is NOT clipped.
%   Table 1 bounds are a marginal geometry screen, not a validated envelope.
%   No heel, leeway, appendage, sail-trim or foil-lift corrections are made.
%   See HULL.md for equations, sources and limitations.

[V, hull, rho, nu, g, geometry] = delft_validate(V, hull, env);
Fn = V ./ sqrt(g * hull.LWL);
moving = V > 0;
% Permit only floating-point roundoff at the two inclusive endpoints.
tol = 8 * eps(0.75);
bad = find(moving & (Fn < 0.15-tol | Fn > 0.75+tol), 1);
if ~isempty(bad)
    error('delft:SpeedRange', ...
        'V(%d) gives Fn=%.16g; require V=0 or 0.15 <= Fn <= 0.75.', ...
        bad, Fn(bad));
end
Re = V .* hull.LWL ./ nu;
bad = find(moving & (~isfinite(Re) | Re <= 100), 1);
if ~isempty(bad)
    error('delft:ReynoldsRange', ...
        'V(%d) gives Re=%.16g; ITTC calculation requires finite Re > 100.', ...
        bad, Re(bad));
end
Cf = nan(size(V));
Rf = zeros(size(V));
Rr = zeros(size(V));
Cf(moving) = 0.075 ./ (log10(Re(moving)) - 2).^2;
Rf(moving) = 0.5 * rho * hull.Swet .* V(moving).^2 .* Cf(moving);
if any(moving(:))
    [nodes, coefficients] = delft_coefficients();
    query = min(max(Fn(moving), nodes(1)), nodes(end));
    a = interp1(nodes, coefficients, query(:), 'linear');
    % Eq. (1.7): a0 is OUTSIDE the volume^(1/3)/LWL multiplier.
    LCB = hull.xFPP - hull.xLCB;
    LCF = hull.xFPP - hull.xLCF;
    terms = [LCB / hull.LWL; hull.Cp; hull.volume^(2/3) / hull.Awp; ...
        hull.BWL / hull.LWL; LCB / LCF; hull.BWL / hull.draft; hull.Cm];
    ratio = a(:,1) + (a(:,2:8) * terms) * hull.volume^(1/3) / hull.LWL;
    Rr(moving) = rho * g * hull.volume .* ratio;
end
R = Rf + Rr;
bad = find(~isfinite(R), 1);
if ~isempty(bad)
    error('delft:NumericalRange', 'Resistance overflow at V(%d). Check SI inputs.', bad);
end
details = struct('Rf', Rf, 'Rr', Rr, 'Fn', Fn, 'Re', Re, 'Cf', Cf);
details.model = 'Keuning-Katgert 2008, eq. 1.7 / table 2; ITTC-1957 with LWL';
details.validity = struct('withinSpeedRange', moving, ...
    'zeroSpeed', ~moving, 'geometry', geometry, ...
    'negativeResiduary', Rr < 0, 'negativeTotal', R < 0, ...
    'reducedHighSpeedDataset', Fn > 0.60, ...
    'note', 'Table 1 marginal screen only; joint and speed-specific validity not established.');
end
