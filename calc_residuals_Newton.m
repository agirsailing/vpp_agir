function [residuals,FM] = calc_residuals_Newton(state, env, moth, TWS, TWA)

% Unpack the 5-DOF state vector requested by the solver
VS         = state(1); % Surge (Forward boat speed)
heave      = state(2); % Vertical displacement
heel       = state(3); % heel [rad]
pitch      = state(4); % pitch [rad]
leeway     = state(5); % drift angle [rad]
s          = state(6); % control 0:2
x_crew     = state(7);   % x-position of crew from origin [m]

[y_crew, alpha_sail] = get_controls(s, moth);

moth.x_cg = (moth.m_boat*moth.x_boat_cg + moth.m_sailor*x_crew) / moth.m_tot;
moth.y_cg = (moth.m_boat*moth.y_boat_cg + moth.m_sailor*y_crew) / moth.m_tot;
moth.r_cg = [moth.x_cg; moth.y_cg; moth.z_cg];

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
[FA, CEA, alfa_sail] = calc_sail(moth, env, alpha_sail, VS, TWS, TWA_eff, heel, pitch);

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
[~, FH_rudder_H, M_rudder_H_local, ~, ~, ~, alfa_elev] = calc_rudder_foil(VS, H_hull, heel, pitch, leeway, moth.rudder_rake, 0, moth.f2_sections, env, moth, false, zeros(6,1));

% Rudder Strut (Vertical)
[~, FH_rudder_V, M_rudder_V_local, ~, ~, ~, alfa_strut] = calc_rudder_foil(VS, H_hull, heel, pitch, leeway, moth.rudder_rake, 0, moth.rv_sections, env, moth, false, zeros(6,1));

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

% --- NEW: HULL SIDE FORCE & INDUCED DRAG ---
% 1. Get the actual submerged draft (T_hull is positive when hull is immersed, 0 when flying)
T_hull = max(0, -H_hull);

% 2. Calculate hull side force using Slender Body Theory (Jones, 1946)
FY_hull = 0.5 * env.water.rho * VS^2 * pi * T_hull^2 * leeway;

% 3. Calculate the induced drag caused by this side force
hull_induced_drag = FY_hull * tan(leeway);

% Convert scalar table outputs into 5-DOF forces and moments
% Add the induced drag penalty to the baseline surge resistance
hull_surge     = -R_hull - hull_induced_drag;            % [N] Resistance acts backwards (-x)
hull_heave     = -env.water.rho * env.g * V_disp;        % [N] Buoyancy acts upwards (-z)
% Buoyancy acts at the hull centreline; its offset from the combined CG gives a roll moment
r_cb = R_boat2earth * ([xB; 0; moth.z_hull_bottom - 0.05] - moth.r_cg);
M_cb = cross(r_cb, [0; 0; hull_heave]);
hull_roll_mom = M_cb(1);

hull_pitch_mom = (xB - moth.x_cg) * (-hull_heave);       % [Nm] Pitch moment (moment arm x buoyancy force)

% Wand Mechanism
% --- WAND MECHANISM (Active Flap Control) ---
% Link the flap angle to the flying height to create a restoring force in heave
dynamic_flap = moth.main_flap_ang + moth.wand_gain * (heave - moth.depth_foil_flying);
dynamic_flap = min(max(dynamic_flap, moth.f1_flap_min), moth.f1_flap_max); % Clamp limits

% Main Foil calculation (Lateral Equilibrium enforcement)
% The vertical foil must produce side force to counter sail, rudder, and hull side force
FY_target = -(FAY + FH_rudder(2) + FY_hull); 
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


% Heel control (s): more s -> more windward -> heel decreases, so g_s rises with s
g_s = (moth.heel_target - heel) / deg2rad(5);
R_s = max(s - 2, min(s - 0, g_s));

% Pitch control (x_crew): moving forward lowers the bow, so g_x rises with x_crew
x_lo = moth.wing_x_pos;  x_hi = moth.wing_x_pos + moth.wing_length;
g_x = (moth.pitch_target - pitch) / deg2rad(2);
R_x = max(x_crew - x_hi, min(x_crew - x_lo, g_x));

residuals = [F_surge; F_heave; M_heel; M_pitch; F_sway; R_s; R_x];


% if nargout > 1
%     AR_s  = moth.b_sail^2/moth.S_sail;
%     CLa_s = 2*pi/(1 + 2*pi/(pi*0.85*AR_s));
%     dg.X               = state(:)';
%     dg.heel_deg        = rad2deg(heel);
%     dg.pitch_deg       = rad2deg(pitch);
%     dg.leeway_deg      = rad2deg(leeway);
%     dg.hull_clearance  = H_hull;
%     dg.sail_CL_uncapped = CLa_s*alfa_sail;       % > CLmax_sail: cap is active
%     dg.FA              = FA';                    % [drive side vertical] N, earth axes
%     dg.FAZ_over_W      = FAZ/moth.weight_tot;    % sail lifting the boat?
%     dg.flap_deg        = rad2deg(dynamic_flap);
%     dg.alpha_main_deg  = rad2deg(out_main.alpha1);
%     dg.main_feasible   = out_main.feasible;
%     dg.alpha_elev_deg  = rad2deg(alfa_elev);
%     dg.alpha_strut_deg = rad2deg(alfa_strut);
%     dg.side_sail_rudder_main = [FAY, FH_rudder(2), F_main(2)];
%     dg.surge_sail_wings_hull_main_rudder = [FAX, F_aero_drag_X, hull_surge, F_main(1), FH_rudder(1)];
%     dg.heel_moment_sail_main_rudder = [M_sail(1), M_main(1), M_rudder(1)];
%     dg.main_drag = out_main.D;
% end

FM.surge_sail = FAX;
FM.surge_wings = F_aero_drag_X;
FM.surge_hull = hull_surge;
FM.surge_main_f = F_main(1);
FM.surge_rudder_f = FH_rudder(1);

FM.heave_sail = FAZ;
FM.heave_hull_float = hull_heave;
FM.heave_main_f = F_main(3);
FM.heave_rudder_f = FH_rudder(3);
FM.heave_hull_weight = moth.weight_tot;

FM.heel_sail = M_sail(1);
FM.heel_hull = hull_roll_mom;
FM.heel_main_f = M_main(1);
FM.heel_rudder_f = M_rudder(1);

FM.pitch_sail = M_sail(2);
FM.pitch_hull = hull_pitch_mom;
FM.pitch_main_f = M_main(2);
FM.pitch_rudder_f = M_rudder(2);
FM.alpha_sail = alpha_sail;
FM.y_crew = y_crew;





% disp("F-surge components: sail="+num2str(FAX)+", wings="+num2str(F_aero_drag_X)+", hull="+num2str(hull_surge)+", main foil="+num2str(F_main(1))+", rudder foil="+num2str(FH_rudder(1))+" [N]")
% disp("F-heave components: sail="+num2str(FAZ)+", hull="+num2str(hull_heave)+", main foil="+num2str(F_main(3))+", rudder foil="+num2str(FH_rudder(2))+", weight="+num2str(moth.weight_tot)+" [N]")
% disp("M_heel components: sail="+num2str(M_sail(1))+", hull="+num2str(hull_roll_mom)+", main foil="+num2str(M_main(1))+", rudder foil="+num2str(M_rudder(1))+" [Nm]")
% disp("M_pitch components: sail="+num2str(M_sail(2))+", hull="+num2str(hull_pitch_mom)+", main foil="+num2str(M_main(2))+", rudder foil="+num2str(M_rudder(2))+" [Nm]")
% disp(" ")

end
