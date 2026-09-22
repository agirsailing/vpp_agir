function [Fbody_N, out] = wings(VairRelBoat_body_mps, moth, env)
%WINGAERORESISTANCE Simple aerodynamic resistance of the Moth wings.
%
% Coordinate system:
%   +x forward, +y starboard, +z downward.
%
% Input:
%   VairRelBoat_body_mps : Air velocity relative to boat [3 x N], m/s
%   moth                 : Structure from get_boat_params()
%   env                  : Structure from get_env_params()
%
% Output:
%   Fbody_N : Aerodynamic force [3 x N], N
%   out     : Useful intermediate results
%
% Model:
%   D = 0.5 * rho * V^2 * CD * S
%
% Source:
%   NASA "Drag Equation"
%   https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/

    % Find the wing entry in moth.windage.
    wingID = find(strcmpi(moth.windage.names, 'wings'), 1);

    % Wing reference area and drag coefficient.
    S  = moth.windage.S(wingID);
    CD = moth.windage.CD(wingID);

    if isnan(S) || isnan(CD)
        error('Set moth.windage.S and CD for the wings.');
    end

    % Air speed.
    speed_mps = sqrt(sum(VairRelBoat_body_mps.^2, 1));

    % Unit vector in the airflow direction.
    flowDirection = zeros(size(VairRelBoat_body_mps));
    moving = speed_mps > 0;

    flowDirection(:, moving) = VairRelBoat_body_mps(:, moving) ./ speed_mps(moving);

    % Dynamic pressure.
    q_Pa = 0.5 .* env.air.rho .* speed_mps.^2;

    % Wing drag magnitude.
    drag_N = q_Pa .* CD .* S;

    % Aerodynamic force vector.
    Fbody_N = flowDirection .* drag_N;

    % Additional outputs.
    out.speed_mps = speed_mps;
    out.dynamicPressure_Pa = q_Pa;
    out.drag_N = drag_N;
    out.longitudinalResistance_N = -Fbody_N(1, :);

end