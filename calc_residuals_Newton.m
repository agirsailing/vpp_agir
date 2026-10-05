function residuals = calc_residuals_Newton(state, env, moth, TW, TWS, TWA)

% Unpack the 5-DOF state vector requested by the solver
VS     = state(1); % Surge (Forward boat speed)
heave  = state(2); % Vertical displacement
heel   = state(3); % Roll angle
pitch  = state(4); % Trim angle
leeway = state(5); % Drift angle

%--------------------------------------------------------------------------
% Aerodynamic calculations
%--------------------------------------------------------------------------
% Adjust true wind angle for lateral drift (leeway)
TWA_eff = TWA - leeway;

% Calculate aerodynamic force components vector FA and Center of Effort
[FA, CEA, ~] = calc_sail(moth, env, TW, VS, TWS, TWA_eff, heel, pitch);

% Divide into force vectors 
FAX = FA(1);
FAY = FA(2);
FAZ = FA(3);

%--------------------------------------------------------------------------
% Hydrodynamic calculations
%--------------------------------------------------------------------------
% Rudder calculation (evaluated first because it uses global leeway)
[~, FH_rudder, M_rudder, ~, ~, ~, ~] = calc_rudder_foil(VS, heave, heel, pitch, leeway, moth.rudder_rake, 0, moth.f2_sections, env, moth, false, zeros(6,1));

% Hull calculation using the pre-calculated table
[R_hull, V_disp, xB] = calc_hull(VS, heave, moth);
% Convert scalar table outputs into 5-DOF forces and moments
hull_surge     = -R_hull;                                % [N] Resistance acts backwards (-x)
hull_heave     = -env.water.rho * env.g * V_disp;        % [N] Buoyancy acts upwards (-z)
hull_roll_mom  = 0;                                      % [Nm] Symmetrical hull has zero roll moment
hull_pitch_mom = (xB - moth.x_cg) * (-hull_heave);       % [Nm] Pitch moment (moment arm x buoyancy force)

% Wand Mechanism
% --- WAND MECHANISM (Active Flap Control) ---
% Link the flap angle to the flying height to create a restoring force in heave
dynamic_flap = moth.main_flap_ang + moth.wand_gain * (heave - moth.depth_foil_flying);
dynamic_flap = min(max(dynamic_flap, moth.f1_flap_min), moth.f1_flap_max); % Clamp limits

% Main Foil calculation (Lateral Equilibrium enforcement)
% The vertical foil must produce side force to counter sail and rudder side force
FY_target = -(FAY + FH_rudder(2)); 
[F_main, M_main, out_main] = calc_main_foil(VS, heave, pitch, dynamic_flap, FY_target, moth, env);


%--------------------------------------------------------------------------
% Equilibrium equations
%--------------------------------------------------------------------------
% Calculate the 5-DOF residuals 
F_surge = FAX + hull_surge + F_main(1) + FH_rudder(1); 
F_heave = FAZ + hull_heave + F_main(3) + FH_rudder(3) + moth.weight_tot; 

% Moments about the boat CG  
M_heel  = FAY * CEA + hull_roll_mom + M_main(1) + M_rudder(1); 
M_pitch = FAX * CEA + hull_pitch_mom + M_main(2) + M_rudder(2);

% Kinematic leeway match 
F_sway  = leeway - out_main.leeway;

% Output array for the Newton solver
residuals = [F_surge; F_heave; M_heel; M_pitch; F_sway];
end
