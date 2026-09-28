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
% M      [Nm] [Mx; My; Mz]  moment about the boat CG -> heel, pitch
%                         (Mx + = heel to stbd, My + = bow up)
% out  truct with diagnostics (angles, CLs, drag, ...)

function [F, M, out] = calc_mainfoil(VS, h, theta, delta, FY, moth, env)

    % Enviroment 
    rho = env.water.rho;                % [kg/m3]
    nu = env.water.nu;                  % [m2/s]
    g = env.g;                          % [m/s2]
    q = 0.5*rho*VS^2;                   % [N/m2] Dynamic pressure
    e = 0.85;                           % [-] Span efficiency factor
    
    % Horizontal foil
    b1 = moth.f1_span;                  % [m] Span
    tc1 = moth.f1_tc;                   % [-] Thickness ratio
    S1 = pi/4*b1*moth.f1_c_root;        % [m2] Planform area of an ellipse
    c1 = S1/b1;                         % [m] Mean chord
    AR1 = b1^2/S1;                      % [-] Aspect ratio
    
    % Vertical
    cs = moth.mv_chord;                 % [m] Chord
    tcs = moth.mv_tc;                   % [-] Thickness ratio
    sw = min(max(h,0), moth.mv_span);   % [m] Wetted vertical span
    Ss = sw*cs;                         % [m2] Wetted vertical area (one side)
    ARs = 2*sw/cs;                      % [-] Effective AR
    
    % Points of application in boat coordinates
    z_foil = moth.z_hull_bottom + moth.mv_span;         % [m] Horizontal foil z
    r_foil = [moth.f1_x - moth.x_cg; 0; z_foil ...
        - moth.z_boat_cg];                              % [m] CG -> foil CE
    r_str = [moth.f1_x - moth.x_cg; 0; z_foil ...
        - sw/2 - moth.z_boat_cg];                       % [m] CG -> vertical CE
    
    % Horizontal foil: lift and drag
    Cla = 2*pi;                                         % [1/rad] 2D lift slope
    alpha0 = moth.f1_alpha0*pi/180;                     % [rad] Zero-lift angle
    alpha1 = theta + delta;                             % [rad] Foil angle of attack
    Cl1 = Cla*(alpha1 - alpha0);                        % [-] 2D lift
    CL1 = Cl1/(1+2/(e*AR1));                            % [-] 3D lift correction
    L1 = q*S1*CL1;                                      % [N] Lift (upwards)
    
    k0 = g/VS^2;                                        % [1/m] Wave number
    CDp1 = 2*ittc(VS*c1/nu)*(1 + 2*tc1 + 60*tc1^4);     % [-] Profile drag (eq. 19-21)
    CDi1 = CL1^2/(pi*e*AR1);                            % [-] Induced drag (eq. 22)
    CDw1 = 0.5*k0*c1*CL1^2*exp(-2*k0*h);                % [-] Wave drag (eq. 28-29)
    D1 = q*S1*(CDp1 + CDi1 + CDw1);                     % [N] Total foil drag
    
    % Vertical foil: side force (sway) and drag
    if sw > 0
        CLs = FY/(q*Ss);                                % [-] Required side force coeff.
        CLas = Cla/(1+2/(e*ARs));                       % [1/rad] 3D lift slope
        leeway = CLs/CLas;                              % [rad] Leeway needed
        CDps = 2*ittc(VS*cs/nu)*(1 + 2*tcs + 60*tcs^4); % [-] Profile drag
        CDis = CLs^2/(pi*e*ARs);                        % [-] Induced drag (eq. 18)
        Dv = q*Ss*(CDps + CDis);                        % [N] Drag
    
        if h < moth.mv_span
            CDspray = 0.009 + 0.013*tcs;                % [-] Spray drag (eq. 27)
            Dspray = q*CDspray*tcs*cs^2;                % [N] 
        else
            Dspray = 0;                                 % Hull in water, no spray
        end
        
    else
        CLs = NaN; leeway = NaN; Dv = 0; Dspray = 0;    % Foil out of the water
    end
    
    % Forces and moments
    F_foil = [-D1; 0; -L1];                             % [N] z down, so lift is negative
    F_str = [-(Dv + Dspray); FY; 0];                    % [N]
    F = F_foil + F_str;                                 % [N] Total force
    M = cross(r_foil,F_foil) + cross(r_str,F_str);      % [Nm] Moment about CG
    
    % Foil output
    out.alpha1 = alpha1;      
    out.CL1 = CL1;
    out.leeway = leeway;      
    out.CLs = CLs;
    out.L1 = L1;
    out.D = struct('foil_profile', q*S1*CDp1, 'foil_induced', q*S1*CDi1, ...
                   'foil_wave', q*S1*CDw1, 'vertical', Dv, 'spray', Dspray);
    out.feasible = sw > 0 && alpha1 >= moth.f1_aoa_min*pi/180 ...
        && alpha1 <= moth.f1_aoa_max*pi/180;
end

function Cf = ittc(Re)
    % ITTC-57 friction line
    Cf = 0.075/(log10(Re)-2)^2;
end
