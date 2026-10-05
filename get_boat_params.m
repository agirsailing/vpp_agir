function moth = get_boat_params()
% origin = stern of hull on deck, x forward, y starboard, z positive DOWN, deck at z = 0
% All moments are taken about the combined CG (boat + sailor): moth.r_cg

    %%  General 
    moth.name        = 'Gungnir';
    
    %%  Mass & crew 
    moth.m_boat      = 55;          % [kg] boat mass
    moth.m_sailor    = 75;          % [kg] sailor
    moth.m_tot       = moth.m_boat+moth.m_sailor;

    moth.x_crew      = 1.500;       % [m] sailor position (fixed LCG)
    moth.z_crew      = NaN;         % [m] sailor CG height (above deck = negative)
    
    %% Center of gravity (boat only)
    moth.x_boat_cg   = 1.770;       % [m]
    moth.y_boat_cg   = 0.00;        % [m]
    moth.z_boat_cg = -0.30;         % [m] vertical CG (above deck)
       
    %%  Hull 
    moth.L_wl        = 3.325;          % [m] waterline length 
    moth.T_hull      = NaN;            % [m] canoe body draught at full displacement
    moth.A_wp        = NaN;            % [m^2] waterplane area
    moth.S_wet_hull  = NaN;            % [m^2] wetted surface at full displacement
    moth.Cp_hull     = NaN;            % [-] prismatic coefficient 
    moth.LCB_hull    = NaN;            % [-] LCB position 
    moth.k_hull      = 0.0;            % [-] hull form factor k, used as (1+k). DSYHS convention: 0
    moth.z_hull_bottom = 0.25;         % [m] deck to hull bottom (z down)
    moth.preCalc = 'preCalcHull.mat';  % [mat] precalculated resistance/displacement table
        
    %%  Main vertical
    moth.mv_profile  = 'NACA0012';
    moth.mv_span     = 1.000;       % [m] hull bottom to main foil 
    moth.mv_chord    = 0.115;       % [m] mean chord
    moth.mv_tc       = 0.12;        % [-] thickness ratio
    
    %%  Main horizontal foil
    moth.f1_profile  = 'NACA2412';
    moth.f1_planform = 'ellipse';   % 'ellipse' or 'trapezoid'
    moth.f1_alpha0   = -2.1;        % [deg] zero-lift angle, from polar
    moth.f1_x        = 1.900;       % [m] CE (centre of effort)
    moth.f1_span     = 1.100;       % [m] 
    moth.f1_c_root   = 0.160;       % [m] max chord
    moth.f1_c_tip    = NaN;         % [m] not used for elliptical planform
    moth.f1_tc       = 0.12;        % [-] thickness ratio
    moth.f1_h_end    = 0;           % [m] end-plate height (0 for none)
    moth.f1_aoa_min  = -3;          % [deg] 
    moth.f1_aoa_max  = 6;           % [deg]
    moth.f1_CLmax    = NaN;         % [-] from polar at aoa_max
    moth.f1_flap_min = deg2rad(-10);% [rad] Physical minimum flap angle
    moth.f1_flap_max = deg2rad(15); % [rad] Physical maximum flap angle
    
    moth.f1_S        = foil_area(moth.f1_planform, moth.f1_span, ...
                                 moth.f1_c_root, moth.f1_c_tip);  % [m^2]
    %% Foil Lift and Drag data
    moth.r_foildata = readtable("DATA_NACA2008.csv");

    %%  Rudder vertical 
    moth.rv_profile  = 'NACA0012';
    moth.rv_span     = 1.100;       % [m] hull level to rudder foil
    moth.rv_chord    = NaN;         % [m] 
    moth.rv_tc       = 0.12;        % [-]
        
    %%  Rudder horizontal foil
    moth.f2_profile  = 'NACA0008';
    moth.f2_planform = 'trapezoid'; % 'ellipse' or 'trapezoid'
    moth.f2_x        = 0.000;       % [m] CE
    moth.f2_span     = 0.600;       % [m]
    moth.f2_c_root   = 0.120;       % [m]
    moth.f2_c_tip    = NaN;         % [m]
    moth.f2_tc       = 0.08;        % [-]
    moth.f2_h_end    = 0;           % [m] 
    moth.f2_aoa_max  = 3;           % [deg] 
    moth.f2_CLmax    = NaN;         % [-]
    moth.f2_sections = [0.00,  0.00, 0.12; 
                        0.00, -0.55, 0.10; 
                        0.00, -1.10, 0.08]'; % [y, z, chord]
    moth.f2_S        = foil_area(moth.f2_planform, moth.f2_span, ...
                                 moth.f2_c_root, moth.f2_c_tip);  % [m^2]
    
    %%  Ride height 
    moth.depth_foil_flying = 0.18;  % [m] main foil depth when foiling
    
    %%  Rig/sail 
    moth.S_sail      = 8.25;        % [m^2] Mach2 sail
    moth.b_sail      = 5.1;         % [m] sail span (AR = b^2/S)
    moth.z_CE_sail   = -2.317;      % [m] sail CE, negative = above hull
    moth.CLmax_sail  = 1.5;         % [-] 
    moth.k_sail      = 1.05;        % [-] sail form factor (1+k)
    
    %% Wings
    moth.numberOfWingbars = 2;       % [-]  
    moth.wingbarSpan_m    = 2.18712; % [m] from the CAD, tip to tip
    moth.wingbarDepth_m   = 0.050;   % [m] from the CAD

    moth.S_wings_m2  = moth.numberOfWingbars*...
        moth.wingbarSpan_m*moth.wingbarDepth_m;
    moth.y_wing_edge = moth.wingbarSpan_m/2; % [m] wing edge from centreline
    moth.wings_CD    = 1.20;

    %%  Derived 
    % Crew lateral position: sitting on the port (windward) wing edge
    moth.y_crew = -moth.y_wing_edge;                                % [m]

    % Combined CG (boat + sailor) = reference point for ALL moments
    moth.x_cg = (moth.m_boat*moth.x_boat_cg + moth.m_sailor*moth.x_crew)/moth.m_tot;  % [m]
    moth.y_cg = (moth.m_boat*moth.y_boat_cg + moth.m_sailor*moth.y_crew)/moth.m_tot;  % [m]
    moth.z_cg = (moth.m_boat*moth.z_boat_cg + moth.m_sailor*moth.z_crew)/moth.m_tot;  % [m]
    moth.r_cg = [moth.x_cg; moth.y_cg; moth.z_cg];                                    % [m]

    moth.x1          = moth.f1_x - moth.x_cg;                       % [m] CG -> foil 1 CE (eq. 17)
    moth.x2          = moth.x_cg - moth.f2_x;                       % [m] CG -> foil 2 CE
    moth.f1_AR       = moth.f1_span^2/moth.f1_S;                    % [-]
    moth.f1_c_mean   = moth.f1_S/moth.f1_span;                      % [m]
    moth.f2_AR       = moth.f2_span^2/moth.f2_S;                    % [-]
    moth.f2_c_mean   = moth.f2_S/moth.f2_span;                      % [m]
    moth.h_fly       = moth.mv_span - moth.depth_foil_flying;       % [m] hull clearance when foiling

end

function S = foil_area(planform, span, c_root, c_tip)
% Planform area of a horizontal foil
    switch lower(planform)
        case 'ellipse'
            S = pi/4*span*c_root;               % [m^2] ellipse, c_root = max chord
        case 'trapezoid'
            S = span*(c_root + c_tip)/2;        % [m^2] straight-tapered
        otherwise
            error('get_boat_params:planform', 'Unknown planform "%s".', planform);
    end
end