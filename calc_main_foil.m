% Forces and moments from the main T-foil (vertical + horizontal foil)
%
% Inputs
% VS     [m/s] boat speed
% h      [m] depth of horizontal foil below free surface (ride height)
% theta  [rad] trim angle, positive = bow up
% delta  [rad] flap angle, modelled as a change of foil incidence
% FY     [N] side force the vertical must produce (sway solved here)
% moth, env - structs from get_boat_params & get_env_params
%
% Outputs (earth-aligned axes: x forward, y starboard, z DOWN)
% F      [N]  [Fx; Fy; Fz]  force on the boat -> surge, sway, heave
% M      [Nm] [Mx; My; Mz]  moment about the COMBINED CG (moth.r_cg)
%                         (Mx + = heel to stbd, My + = bow up)
% out    struct with diagnostics (angles, CLs, drag, ...)

function [F, M, out] = calc_main_foil(VS, h, heel, theta, delta, FY, moth, env)

    % Enviroment 
    rho = env.water.rho;                % [kg/m3]
    g = env.g;                          % [m/s2]
    q = 0.5*rho*VS^2;                   % [N/m2] Dynamic pressure
    e = 0.85;                           % [-] Span efficiency factor
    Cla = 2*pi;                         % [1/rad] 2D lift slope
    
    % Horizontal foil (geometry from get_boat_params)
    S1 = moth.f1_S;                     % [m2] Planform area
    c1 = moth.f1_c_mean;                % [m] Mean chord
    AR1 = moth.f1_AR;                   % [-] Aspect ratio
    
    % Vertical
    cs = moth.mv_chord;                 % [m] Chord
    tcs = moth.mv_tc;                   % [-] Thickness ratio
    sw = min(max(h,0), moth.mv_span);   % [m] Wetted vertical span
    Ss = sw*cs;                         % [m2] Wetted vertical area (one side)
    ARs = 2*sw/cs;                      % [-] Effective AR
    
    % % Points of application, relative to the combined CG
    % z_foil = moth.z_hull_bottom + moth.mv_span;         % [m] Horizontal foil z
    % r_foil = [moth.f1_x; 0; z_foil] - moth.r_cg;        % [m] CG -> foil CE
    % r_str = [moth.f1_x; 0; z_foil - sw/2] - moth.r_cg;  % [m] CG -> vertical CE
    
    % Points of application, relative to the combined CG (Rotated to Earth Frame)
    z_foil = moth.z_hull_bottom + moth.mv_span;

    % Use 'heel' and the existing 'theta' (pitch)
    R_heel  = [1 0 0; 0 cos(heel) -sin(heel); 0 sin(heel) cos(heel)];
    R_pitch = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    R_boat2earth = R_pitch * R_heel;
    
    r_foil = R_boat2earth * ([moth.f1_x; 0; z_foil] - moth.r_cg);
    r_str = R_boat2earth * ([moth.f1_x; 0; z_foil - sw/2] - moth.r_cg);

    % Horizontal foil: lift and drag
    alpha1 = theta + delta;                             % [rad] Foil angle of attack
    foil_wet = h > 0;                                   % Foil below the free surface
    
    if foil_wet
        alpha0 = moth.f1_alpha0*pi/180;                 % [rad] Zero-lift angle
        Cl1 = Cla*(alpha1 - alpha0);                    % [-] 2D lift
        CL1 = Cl1/(1+2/(e*AR1));                        % [-] 3D lift correction
        L1 = q*S1*CL1;                                  % [N] Lift (upwards)
        
        k0 = g/VS^2;                                    % [1/m] Wave number
        CDp1 = 0.008 + 0.01*Cl1^2;                      % [-] Section drag (course VPP)
        CDi1 = CL1^2/(pi*e*AR1);                        % [-] Induced drag (eq. 22)
        CDw1 = 0.5*k0*c1*CL1^2*exp(-2*k0*h);            % [-] Wave drag (eq. 28-29)
    else
        CL1 = 0; L1 = 0;                                % Foil out of the water
        CDp1 = 0; CDi1 = 0; CDw1 = 0;
    end
    D1 = q*S1*(CDp1 + CDi1 + CDw1);                     % [N] Total foil drag
    
    % Earth-y force of the foil = cos(heel)*FYloc + sin(heel)*L1 must equal FY,
% so the strut only supplies what the tilted lift does not
FYloc = (FY - sin(heel)*L1)/cos(heel);

    % Vertical foil: side force (sway) and drag
    if sw > 0
        % CLs = FY/(q*Ss); 
        CLs = FYloc/(q*Ss);                               % [-] Required side force coeff. (3D)
        CLas = Cla/(1+2/(e*ARs));                       % [1/rad] 3D lift slope
        leeway = CLs/CLas;                              % [rad] Leeway needed
        Cls = Cla*leeway;                               % [-] 2D section lift
        CDps = 0.008 + 0.01*Cls^2;                      % [-] Section drag (course VPP)
        CDis = CLs^2/(pi*e*ARs);                        % [-] Induced drag (eq. 18)
        Dv = q*Ss*(CDps + CDis);                        % [N] Drag
    
        if h < moth.mv_span
            CDspray = 0.009 + 0.013*tcs;                % [-] Spray drag (eq. 27)
            Dspray = q*CDspray*tcs*cs^2;                % [N] 
        else
            Dspray = 0;                                 % Hull in water, no spray
        end
        
    else
        CLs = NaN; leeway = NaN; Dv = 0; Dspray = 0;    % Vertical out of the water
    end
    
    % Forces and moments
    F_foil = R_boat2earth*[-D1; 0; -L1];                             % [N] z down, so lift is negative
    F_str  = R_boat2earth*[-(Dv + Dspray); FYloc*(sw > 0); 0];           % [N] no side force if dry
    F = F_foil + F_str;                                 % [N] Total force
    M = cross(r_foil,F_foil) + cross(r_str,F_str);      % [Nm] Moment about combined CG
    
    % Foil output
    out.alpha1 = alpha1;      
    out.CL1 = CL1;
    out.leeway = leeway;      
    out.CLs = CLs;
    out.L1 = L1;
    out.D = struct('foil_profile', q*S1*CDp1, 'foil_induced', q*S1*CDi1, ...
                   'foil_wave', q*S1*CDw1, 'vertical', Dv, 'spray', Dspray);
    out.feasible = foil_wet && alpha1 >= moth.f1_aoa_min*pi/180 ...
        && alpha1 <= moth.f1_aoa_max*pi/180;
end