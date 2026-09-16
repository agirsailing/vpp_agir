function moth = get_boat_params()
% origin = stern of hull on deck, x forward, z positive DOWN, deck at z = 0

    %%  General 
    moth.name        = 'Gungnir';
    
    %%  Mass & crew 
    moth.m_boat      = 50;          % [kg] boat mass
    moth.m_sailor    = 75;          % [kg] sailor
    moth.x_boat_cg   = 1.770;       % [m] boat CG from origin
    moth.x_crew      = 1.500;       % [m] sailor position (fixed LCG)
    moth.y_wing_edge = 1.125;       % [m] wing edge from centreline
   
    %%  Hull 
    moth.L_wl        = 3.325;       % [m] waterline length (LOA used until measured)
    moth.B_wl        = NaN;         % [m] waterline beam 
    moth.T_hull      = NaN;         % [m] canoe body draught at full displacement
    moth.A_wp        = NaN;         % [m^2] waterplane area
    moth.S_wet_hull  = NaN;         % [m^2] wetted surface at full displacement
    moth.Cp_hull     = NaN;         % [-] prismatic coefficient 
    moth.LCB_hull    = NaN;         % [-] LCB position 
    moth.k_hull      = NaN;         % [-] hull form factor (1+k)
    
    %%  Main vertical
    moth.mv_profile  = 'NACA0012';
    moth.mv_span     = 1.000;       % [m] hull bottom to main foil 
    moth.mv_chord    = NaN;         % [m] mean chord
    moth.mv_tc       = 0.12;        % [-] thickness ratio
    
    %%  Main horizontal foil
    moth.f1_profile  = 'NACA2412';
    moth.f1_x        = 1.900;       % [m] CE (centre of effort)
    moth.f1_span     = 1.100;       % [m] 
    moth.f1_c_root   = 0.160;       % [m] max chord
    moth.f1_c_tip    = NaN;         % [m] elliptical planform
    moth.f1_tc       = 0.12;        % [-] thickness ratio
    moth.f1_h_end    = 0;           % [m] end-plate height (0 for none)
    moth.f1_aoa_min  = -3;          % [deg] 
    moth.f1_aoa_max  = 6;           % [deg]
    moth.f1_CLmax    = NaN;         % [-] from polar at aoa_max
    
    %%  Rudder vertical 
    moth.rv_profile  = 'NACA0012';
    moth.rv_span     = 1.100;       % [m] hull level to rudder foil
    moth.rv_chord    = NaN;         % [m] 
    moth.rv_tc       = 0.12;        % [-]
        
    %%  Rudder horizontal foil
    moth.f2_profile  = 'NACA0008';
    moth.f2_x        = 0.000;       % [m] CE
    moth.f2_span     = 0.600;       % [m]
    moth.f2_c_root   = NaN;         % [m]
    moth.f2_c_tip    = NaN;         % [m]
    moth.f2_S        = NaN;         % [m^2]
    moth.f2_tc       = 0.08;        % [-]
    moth.f2_h_end    = 0;           % [m] 
    moth.f2_aoa_max  = 3;           % [deg] 
    moth.f2_CLmax    = NaN;         % [-]
    
    %%  Ride height 
    moth.depth_foil_flying = 0.18;  % [m] main foil depth when foiling
    
    %%  Rig/sail 
    moth.S_sail      = 8.25;        % [m^2] Mach2 sail
    moth.b_sail      = NaN;         % [m] sail span (AR = b^2/S)
    moth.z_CE_sail   = -1.690;      % [m] sail CE, negative = above hull
    moth.CLmax_sail  = 1.5;         % [-] 
    moth.k_sail      = 1.05;        % [-] sail form factor (1+k)
    
    %%  Windage 
    moth.windage.names = {'hull','wings','crew','rig','gantry'};
    moth.windage.S     = [NaN NaN NaN NaN NaN];                     % [m^2]
    moth.windage.CD    = [NaN NaN NaN NaN NaN];                     % [-]
    
    %%  Derived 
    moth.m_total     = moth.m_boat + moth.m_sailor;                 % [kg]
    moth.x_cg        = (moth.m_boat*moth.x_boat_cg + ...
        moth.m_sailor*moth.x_crew)/moth.m_total;                    % [m]
    moth.x_RM        = moth.y_wing_edge + moth.y_hike;              % [m] righting arm (eq. 2)
    moth.x1          = moth.f1_x - moth.x_cg;                       % [m] CG -> foil 1 CE (eq. 17)
    moth.x2          = moth.x_cg - moth.f2_x;                       % [m] CG -> foil 2 CE
    moth.f1_AR       = moth.f1_span^2/moth.f1_S;                    % [-] (8.66)
    moth.f1_c_mean   = moth.f1_S/moth.f1_span;                      % [m]
    moth.f2_AR       = moth.f2_span^2/moth.f2_S;                    % [-]
    moth.f2_c_mean   = moth.f2_S/moth.f2_span;                      % [m]
    moth.h_fly       = moth.mv_span - moth.depth_foil_flying;       % [m] hull clearance when foiling

end