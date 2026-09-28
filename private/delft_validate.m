function [V, h, rho, nu, g, geometry] = delft_validate(V, h, env)
%DELFT_VALIDATE Input contract and marginal ranges of 2008 Table 1.
if ~(isnumeric(V) && isreal(V) && isvector(V) && ~isempty(V))
    error('delft:InvalidSpeed', 'V must be a nonempty real numeric scalar or vector.');
end
V = double(V);
bad = find(~isfinite(V) | V < 0, 1);
if ~isempty(bad)
    error('delft:InvalidSpeed', 'V(%d) must be finite and nonnegative.', bad);
end
checkStruct(h, 'hull');
positive = {'LWL','BWL','draft','volume','Swet','Awp','Cp','Cm'};
coordinates = {'xLCB','xLCF','xFPP'};
for field = [positive coordinates]
    name = field{1};
    h.(name) = scalarField(h, name, ['hull.' name]);
    if ismember(name, positive) && h.(name) <= 0
        error('delft:InvalidInput', 'hull.%s must be positive.', name);
    end
end
checkStruct(env, 'env');
if ~isfield(env, 'water')
    error('delft:MissingField', 'Missing env.water.');
end
checkStruct(env.water, 'env.water');
rho = scalarField(env.water, 'rho', 'env.water.rho');
nu = scalarField(env.water, 'nu', 'env.water.nu');
g = scalarField(env, 'g', 'env.g');
values = [rho nu g];
names = {'env.water.rho','env.water.nu','env.g'};
for k = 1:3
    if values(k) <= 0
        error('delft:InvalidInput', '%s must be positive.', names{k});
    end
end
LCB = h.xFPP - h.xLCB;
LCF = h.xFPP - h.xLCF;
if ~(LCB > 0 && LCB < h.LWL && LCF > 0 && LCF < h.LWL)
    error('delft:InvalidGeometry', ...
        'hull.xLCB and hull.xLCF must lie between hull.xFPP-hull.LWL and hull.xFPP.');
end
% Mesh integration can put exact physical boundaries a few ulps above one.
physicalRatios = [h.Cp h.Cm h.Awp/h.LWL/h.BWL ...
    h.volume/h.LWL/h.BWL/h.draft];
if any(physicalRatios > 1 + 8*eps(max(1,abs(physicalRatios))))
    error('delft:InvalidGeometry', ...
        'Check hull.Cp, hull.Cm, hull.Awp and hull.volume against hull dimensions.');
end
% Extrema of printed (rounded) columns in Keuning-Katgert 2008 Table 1.
% Independent bounds do not describe the joint domain or speed dependence.
labels = {'LCB_fpp/LWL','Cp','volume^(2/3)/Awp','BWL/LWL', ...
    'LCB_fpp/LCF_fpp','volume^(1/3)/LWL','Cm','BWL/draft'};
ratios = [LCB/h.LWL, h.Cp, h.volume^(2/3)/h.Awp, h.BWL/h.LWL, ...
    LCB/LCF, h.volume^(1/3)/h.LWL, h.Cm, h.BWL/h.draft];
lower = [0.500 0.519 0.079 0.170 0.920 0.12 0.646 2.46];
upper = [0.582 0.599 0.265 0.366 1.002 0.23 0.790 19.38];
tol = 8 * eps(max(1, abs(ratios)));
bad = find(~isfinite(ratios) | ratios < lower-tol | ratios > upper+tol, 1);
if ~isempty(bad)
    error('delft:GeometryRange', ...
        'hull ratio %s=%.16g is outside the Table 1 screen [%.6g, %.6g].', ...
        labels{bad}, ratios(bad), lower(bad), upper(bad));
end
geometry = struct('passed', true, 'basis', 'Keuning-Katgert 2008 Table 1 printed marginal extrema', ...
    'names', {labels}, 'values', ratios, 'lower', lower, 'upper', upper);
end

function checkStruct(value, label)
if ~(isstruct(value) && isscalar(value))
    error('delft:InvalidInput', '%s must be a scalar structure.', label);
end
end

function value = scalarField(s, field, label)
if ~isfield(s, field)
    error('delft:MissingField', 'Missing %s.', label);
end
value = s.(field);
if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
    error('delft:InvalidInput', '%s must be a finite real numeric scalar.', label);
end
value = double(value);
end
