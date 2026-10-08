function residuals = calc_residuals_Newton(state, env, moth, TW, TWS, TWA)

% Unpack the 5-DOF state vector requested by the solver
VS     = state(1); % Surge (Forward boat speed)
heave  = state(2); % Vertical displacement
heel   = state(3); % Roll angle
pitch  = state(4); % Trim angle
leeway = state(5); % Drift angle

H_hull = moth.mv_span - heave; % Negative when hull is immersed

R_heel  = [1 0 0; 0 cos(heel) -sin(heel); 0 sin(heel) cos(heel)];
R_pitch = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
R_boat2earth = R_pitch * R_heel;

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
% Parasitic Aerodynamic Drag (Wings + Windage)
%--------------------------------------------------------------------------
% Apparent wind velocity vector in boat frame: [Vx; Vy; Vz]
V_air_rel = [ -VS - TWS*cos(TWA_eff);
              -TWS*sin(TWA_eff);
               0 ];

% Wing aero drag via calc_wings.m
[F_wings, ~] = calc_wings(V_air_rel, moth, env);

% Total aero resistance in surge (-x direction)
F_aero_drag_X = F_wings(1);

%--------------------------------------------------------------------------
% Hydrodynamic calculations
%--------------------------------------------------------------------------
% Rudder Elevator (Horizontal)
[~, FH_rudder_H, M_rudder_H_local, ~, ~, ~, ~] = calc_rudder_foil(VS, H_hull, heel, pitch, leeway, moth.rudder_rake, 0, moth.f2_sections, env, moth, false, zeros(6,1));

% Rudder Strut (Vertical)
[~, FH_rudder_V, M_rudder_V_local, ~, ~, ~, ~] = calc_rudder_foil(VS, H_hull, heel, pitch, leeway, moth.rudder_rake, 0, moth.rv_sections, env, moth, false, zeros(6,1));

% Combine local forces and moments
FH_rudder = FH_rudder_H + FH_rudder_V;
M_rudder_local = M_rudder_H_local + M_rudder_V_local;

% Transfer total rudder moments to combined CG (moth.r_cg)
% The local origin (0,0,0) of the rudder panel coordinates is the hull attachment point
r_rudder = R_boat2earth * ([moth.f2_x; 0; moth.z_hull_bottom] - moth.r_cg);
M_rudder = cross(r_rudder, FH_rudder) + M_rudder_local;

% Sail moment about combined CG
r_sail = R_boat2earth * ([moth.x_boat_cg; 0; moth.z_CE_sail] - moth.r_cg);
M_sail = cross(r_sail, FA);

% Hull calculation using the pre-calculated table
[R_hull, V_disp, xB] = calc_hull(VS, H_hull, moth);

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
[F_main, M_main, out_main] = calc_main_foil(VS, heave, heel, pitch, dynamic_flap, FY_target, moth, env);


%--------------------------------------------------------------------------
% Equilibrium equations
%--------------------------------------------------------------------------
% Calculate the 5-DOF residuals 
F_surge = FAX + F_aero_drag_X +  hull_surge + F_main(1) + FH_rudder(1); 
F_heave = FAZ + hull_heave + F_main(3) + FH_rudder(3) + moth.weight_tot; 

% Moments about the boat CG  
M_heel  = M_sail(1) + hull_roll_mom + M_main(1) + M_rudder(1); 
M_pitch = M_sail(2) + hull_pitch_mom + M_main(2) + M_rudder(2);

% Kinematic leeway match 
F_sway  = leeway - out_main.leeway;


% Output array for the Newton solver
residuals = [F_surge; F_heave; M_heel; M_pitch; F_sway];
end
