%----------------------------------------------------------------  
% Initialize
%----------------------------------------------------------------  
clear
clear global
close all

addpath("JavaFoil\")
addpath("Foils_data\")
addpath("Helper_functions\")

global DATA % A test-variable for Razola & Kuttenkeuler only :-)
 
%---------------------------------------------------------------  
% Indata
%----------------------------------------------------------------  

% Load parameters for the moth and environment
env = get_env_params();
moth = get_boat_params();
disp(moth.r_cg')
assert(all(isfinite(moth.r_cg)),'r_cg contains NaN');

% total weight of moth
moth.weight_tot = moth.m_tot * env.g;

% VPP control settings
moth.main_flap_ang = 0.0; % [rad] Flap angle on horizontal main foil
moth.wand_gain = 0.8;    % [rad/m] "Gearing": Flap sensitivity to heave

moth.rudder_rake   = deg2rad(0.0); % [rad] Rudder rake angle



%--------------------------------------------------------------------------

H = -0.5; % [m] Constant flying height

% Define vectors for the required homework span
TWS_vec = 2.5:1:7.5;      % [m/s] True wind speeds
TWA_vec = 5:5:180;         % [deg] True wind angles (10-degree steps)

nTWS = length(TWS_vec);
nTWA = length(TWA_vec);

% Pre-allocate results matrices
VS_res     = NaN(nTWS, nTWA);   % surge speed [m/s]
heave_res  = NaN(nTWS, nTWA);   % heave [m]
heel_res   = NaN(nTWS, nTWA);   % heel [deg]
pitch_res  = NaN(nTWS, nTWA);   % pitch [deg]
leeway_res = NaN(nTWS, nTWA);   % leeway [deg]
alpha_res  = NaN(nTWS, nTWA);   % sail angle
y_crew_res = NaN(nTWS, nTWA);   % y crew position [m]
x_crew_res = NaN(nTWS, nTWA);   % x crew position [m]

% Surge forces [N]
surge_sail_res    = NaN(nTWS, nTWA);
surge_wings_res   = NaN(nTWS, nTWA);
surge_hull_res    = NaN(nTWS, nTWA);
surge_main_res    = NaN(nTWS, nTWA);
surge_rudder_res  = NaN(nTWS, nTWA);

% Heave forces [N]
heave_sail_res    = NaN(nTWS, nTWA);
heave_hull_res    = NaN(nTWS, nTWA);
heave_main_res    = NaN(nTWS, nTWA);
heave_rudder_res  = NaN(nTWS, nTWA);
heave_weight_res  = NaN(nTWS, nTWA);

% Heel moments [Nm]
heel_sail_res     = NaN(nTWS, nTWA);
heel_main_res     = NaN(nTWS, nTWA);
heel_rudder_res   = NaN(nTWS, nTWA);
heel_crew_res     = NaN(nTWS, nTWA);
heel_hull_res     = NaN(nTWS, nTWA);

% Pitch moments [Nm]
pitch_sail_res    = NaN(nTWS, nTWA);
pitch_hull_res    = NaN(nTWS, nTWA);
pitch_main_res    = NaN(nTWS, nTWA);
pitch_rudder_res  = NaN(nTWS, nTWA);
pitch_crew_res    = NaN(nTWS, nTWA);

% sail aero
sail_CL_res       = NaN(nTWS, nTWA);
sail_CD_res       = NaN(nTWS, nTWA);
sail_awa_res      = NaN(nTWS, nTWA);

% Run the sweeps
for i = 1:nTWS
    TWS = TWS_vec(i);
    fprintf('\n====================================\n');
    fprintf('Starting Wind Speed: %d m/s\n', TWS);

    moth.heel_target = moth.heel_target_foiling;
    mh = moth;  mh.heel_target = moth.heel_target_disp;       % hull-borne settings
    mf = moth;  mf.heel_target = moth.heel_target_foiling;    % foiling settings
    Xh0 = [max(TWS*0.5,1.5); moth.mv_span+0.10; 0; 0; deg2rad(-2); 1; 1.5];
    Xf0 = [9; moth.depth_foil_flying; deg2rad(-20); 0; deg2rad(-0.3); 1; 1.5];
    
    % % 5-DOF initial guess: [VS; heave; heel; pitch; leeway]
    % X0 = [TWS * 1.3; 0.88; deg2rad(0); deg2rad(1.0); deg2rad(-3.0)];
    %
    %
    for j = 1:nTWA
        TWA_deg = TWA_vec(j);
        TWA_rad = deg2rad(TWA_deg); 
        
        fprintf('  TWA %3d deg... ', TWA_deg);

        %  Dynamic Heel Target 
        % (cos(TWA_rad) + 1) / 2 creates a smooth multiplier:
        % TWA 0   -> Multiplier = 1.0 (Full Veal heel)
        % TWA 90  -> Multiplier = 0.5 (Half Veal heel)
        % TWA 180 -> Multiplier = 0.0 (Flat boat)
        heel_fade = (cos(TWA_rad) + 1) / 2; 
        
        mh.heel_target = moth.heel_target_disp * heel_fade;
        mf.heel_target = moth.heel_target_foiling * heel_fade;

       % Call the 5DOF Newton solver
        % 1) Hull-borne solution
        [Xh, itH, FLAGh, FMh] = solve_Newton(Xh0, env, mh, TWS, TWA_rad);
        disp_ok = (FLAGh == 0 && Xh(1) >= 0.5 && Xh(2) > moth.mv_span);   
        
        % Always accept a converged displacement solution as a safe baseline
        if FLAGh == 0
            X = Xh; iter = itH; FLAG = FLAGh; FM = FMh;
            Xh0 = Xh; 
        else
            X = Xh; FLAG = FLAGh; % pass failure state forward
        end

        % 2) Attempt Foiling if hull lifted out OR displacement failed
        if ~disp_ok || FLAGh ~= 0
            [Xf, itF, FLAGf, FMf] = solve_Newton(Xf0, env, mf, TWS, TWA_rad);
            if FLAGf == 0 && Xf(1) >= 1
                X = Xf; iter = itF; FLAG = FLAGf; FM = FMf;
                Xf0 = Xf;                              
            elseif FLAGh == 0
                % Foiling failed, but we have a valid displacement baseline. 
                % Revert to the Xh solution and DO NOT flag as a failure.
                X = Xh; iter = itH; FLAG = FLAGh; FM = FMh;
            else
                FLAG = 4; % Both solvers completely failed
            end
        end

        VS = X(1);
        heave = X(2);
        heel = rad2deg(X(3));
        pitch = rad2deg(X(4));
        leeway = rad2deg(X(5));
        [y_crew, alpha_now] = get_controls(X(6), moth);
        x_crew = X(7);


        % Safety Check: Did the point fail, stall, or hit the 50 iteration limit?
        if FLAG ~= 0 || isnan(VS) || VS < 1 || (moth.mv_span - X(2)) < -0.16 || abs(X(5)) > deg2rad(8)
            fprintf('Stalled/Failed (FLAG %d). Resetting.\n', FLAG);
            disp("Surge="+num2str(VS)+"m/s, Heave="+num2str(heave)+"m, Heel="+num2str(heel)+"deg, Pitch="+num2str(pitch)+"deg, leeway="+num2str(leeway)+"deg, s="+num2str(X(6))+".")
            VS_res(i,j) = NaN; % Record as 0 to indicate it cannot foil here
            heave_res(i,j)     = NaN;
            heel_res(i,j)      = NaN;
            pitch_res(i,j)     = NaN;
            leeway_res(i,j)    = NaN;
            alpha_res(i,j)     = NaN;
            y_crew_res(i,j)    = NaN;
            x_crew_res(i,j)    = NaN;


        else
            fprintf('Success! VS = %5.2f m/s (iters: %d)\n', VS, iter);
            VS_res(i,j) = VS;

            heave_res(i,j)  = heave;
            heel_res(i,j)   = heel;
            pitch_res(i,j)  = pitch;
            leeway_res(i,j) = leeway;
            alpha_res(i,j) = rad2deg(alpha_now);
            y_crew_res(i,j) = y_crew;
            x_crew_res(i,j) = x_crew;
            % Aerodynamic sail coefficients
            sail_CL_res(i,j)  = FM.sail_CL;
            sail_CD_res(i,j)  = FM.sail_CD;
            sail_awa_res(i,j) = FM.AppWinAng;
            % Surge components
            surge_sail_res(i,j)   = FM.surge_sail;
            surge_wings_res(i,j)  = FM.surge_wings;
            surge_hull_res(i,j)   = FM.surge_hull;
            surge_main_res(i,j)   = FM.surge_main_f;
            surge_rudder_res(i,j) = FM.surge_rudder_f;

            % Heave components
            heave_sail_res(i,j)   = FM.heave_sail;
            heave_hull_res(i,j)   = FM.heave_hull_float;
            heave_main_res(i,j)   = FM.heave_main_f;
            heave_rudder_res(i,j) = FM.heave_rudder_f;
            heave_weight_res(i,j) = FM.heave_hull_weight;

            % Heel components
            heel_sail_res(i,j)    = FM.heel_sail;
            heel_main_res(i,j)    = FM.heel_main_f;
            heel_rudder_res(i,j)  = FM.heel_rudder_f;
            heel_crew_res(i,j)    = FM.heel_crew;
            heel_hull_res(i,j)    = FM.heel_hull;

            % Pitch components
            pitch_sail_res(i,j)   = FM.pitch_sail;
            pitch_hull_res(i,j)   = FM.pitch_hull;
            pitch_main_res(i,j)   = FM.pitch_main_f;
            pitch_rudder_res(i,j) = FM.pitch_rudder_f;
            pitch_crew_res(i,j)   = FM.pitch_crew;

            disp("Surge="+num2str(VS)+"m/s, Heave="+num2str(heave)+"m, Heel="+num2str(heel)+"deg, Pitch="+num2str(pitch)+"deg, leeway="+num2str(leeway)+"deg, s="+num2str(X(6))+".")
            % Update guess for the next angle (Sequential Updating)
        end
    end
end


%%
%--------------------------------------------------------------------------
% Polar plot (boat speed vs TWA, one curve per TWS)
%--------------------------------------------------------------------------
cmap = parula(nTWS + 2);
cmap = cmap(1:nTWS, :);

figure('Name','Polar','Color','w','Position',[100 100 800 650]);
for k = 1:nTWS
    polarplot(deg2rad(TWA_vec), VS_res(k,:), '-', 'Color', cmap(k,:), 'LineWidth', 2, ...
              'Marker','o', 'MarkerSize', 4, 'MarkerFaceColor', cmap(k,:), ...
              'DisplayName', sprintf('TWS = %g m/s', TWS_vec(k)));
    hold on
end
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.ThetaTick = 0:30:180;
pax.RAxisLocation = 0;
pax.GridAlpha = 0.2;
pax.FontSize = 10;
thetalim([0 180]);
rlim([0 max(15, ceil(max(VS_res(:),[],'omitnan')))]);
title(sprintf('Velocity prediction polar'),'FontWeight','bold');
lg = legend('Location','eastoutside'); lg.Box = 'off';

% Rotated radial-axis title, left of the tick labels (placed last)
drawnow;
ax_pos = pax.Position;
% Rotated radial-axis title, left of the tick labels on the 0 deg line
rmax = pax.RLim(2);
text(pax, deg2rad(-9), 0.5*rmax, 'Boat speed [m/s]', ...
     'Rotation', 90, 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', 'FontSize', 11, 'Clipping', 'off');
%%
%--------------------------------------------------------------------------
% State variables vs TWA (one curve per TWS, same colours as the polar)
%--------------------------------------------------------------------------
panels = { VS_res,     'Boat speed [m/s]';
           heave_res,  'Heave [m]';
           heel_res,   'Heel [deg]';
           pitch_res,  'Pitch [deg]';
           leeway_res, 'Leeway [deg]';
           alpha_res,  'Sail AoA [deg]' };

fig = figure('Name','States','Color','w','Position',[100 100 1200 850]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Solved state vs true wind angle', 'FontWeight','bold');

for p = 1:size(panels,1)
    ax = nexttile; hold(ax,'on');
    for k = 1:nTWS
        plot(ax, TWA_vec, panels{p,1}(k,:), '-', 'Color', cmap(k,:), 'LineWidth', 1.8, ...
             'DisplayName', sprintf('TWS = %g m/s', TWS_vec(k)));
    end
    grid(ax,'on');
    set(ax, 'Box','on', 'FontSize',10, 'GridAlpha',0.15, 'LineWidth',0.8, ...
            'XLim',[min(TWA_vec) max(TWA_vec)]);
    ylabel(ax, panels{p,2});
    if p > 6, xlabel(ax, 'TWA [deg]'); end
    if p == 4; ylim(ax, [-1 1]);end
end

% one shared legend on the right
lg = legend(ax, 'Location','eastoutside'); lg.Layout.Tile = 'east'; lg.Box = 'off';
%%
%--------------------------------------------------------------------------
% Force and moment balances (one figure per TWS)
%--------------------------------------------------------------------------
col.sail   = [0.00 0.45 0.74];
col.main   = [0.85 0.33 0.10];
col.rudder = [0.49 0.18 0.56];
col.hull   = [0.30 0.75 0.93];
col.crew   = [0.47 0.67 0.19];
col.wings  = [0.93 0.69 0.13];
col.weight = [0.40 0.40 0.40];
col.sum    = [0 0 0];
lw = 1.8;

% helper: plot one series with consistent style
pl = @(y, c, name) plot(TWA_vec, y, '-', 'Color', c, 'LineWidth', lw, 'DisplayName', name);

% helper: common axes formatting
fmt = @(ttl, yl) set_axes(ttl, yl, TWA_vec);

for i = 1:nTWS
    curr_tws_kn = TWS_vec(i) * env.ms2kn;

    fig = figure('Name', sprintf('Force & Moment Breakdown - TWS %.1f kn', curr_tws_kn), ...
                 'Color', 'w', 'Position', [100 100 1200 800]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Force and moment balance, TWS = %g m/s (%.1f kn)', TWS_vec(i), curr_tws_kn), ...
          'FontWeight', 'bold');

    % 1. Surge
    nexttile; hold on;
    pl(surge_sail_res(i,:),   col.sail,   'Sail drive');
    pl(surge_wings_res(i,:),  col.wings,  'Wings drag');
    pl(surge_main_res(i,:),   col.main,   'Main foil drag');
    pl(surge_rudder_res(i,:), col.rudder, 'Rudder drag');
    pl(surge_hull_res(i,:),   col.hull,   'Hull drag');
    s_surge = surge_sail_res(i,:) + surge_wings_res(i,:) + surge_main_res(i,:) ...
            + surge_rudder_res(i,:) + surge_hull_res(i,:);
    plot(TWA_vec, s_surge, '--', 'Color', col.sum, 'LineWidth', 1.4, 'DisplayName', 'Sum (should be ~0)');
    fmt('Surge forces (X)', 'Force [N]');

    % 2. Heave
    nexttile; hold on;
    pl(heave_sail_res(i,:),   col.sail,   'Sail vertical force');
    pl(heave_main_res(i,:),   col.main,   'Main foil lift');
    pl(heave_rudder_res(i,:), col.rudder, 'Rudder elevator lift');
    pl(heave_hull_res(i,:),   col.hull,   'Hull buoyancy');
    pl(heave_weight_res(i,:), col.weight, 'Total weight');
    s_heave = heave_sail_res(i,:) + heave_main_res(i,:) + heave_rudder_res(i,:) ...
            + heave_hull_res(i,:) + heave_weight_res(i,:);
    plot(TWA_vec, s_heave, '--', 'Color', col.sum, 'LineWidth', 1.4, 'DisplayName', 'Sum (should be ~0)');
    fmt('Heave forces (Z, + = down)', 'Force [N]');

    % 3. Heel moment
    nexttile; hold on;
    pl(heel_sail_res(i,:),   col.sail,   'Sail');
    pl(heel_main_res(i,:),   col.main,   'Main foil');
    pl(heel_rudder_res(i,:), col.rudder, 'Rudder');
    pl(heel_hull_res(i,:),   col.hull,   'Hull buoyancy');
    pl(heel_crew_res(i,:),   col.crew,   'Crew weight');
    s_heel = heel_sail_res(i,:) + heel_main_res(i,:) + heel_rudder_res(i,:) ...
           + heel_hull_res(i,:) + heel_crew_res(i,:);
    plot(TWA_vec, s_heel, '--', 'Color', col.sum, 'LineWidth', 1.4, 'DisplayName', 'Sum (should be ~0)');
    fmt('Heeling moments (Mx about boat CG)', 'Moment [Nm]');

    % 4. Pitch moment
    nexttile; hold on;
    pl(pitch_sail_res(i,:),   col.sail,   'Sail');
    pl(pitch_main_res(i,:),   col.main,   'Main foil');
    pl(pitch_rudder_res(i,:), col.rudder, 'Rudder');
    pl(pitch_hull_res(i,:),   col.hull,   'Hull buoyancy');
    pl(pitch_crew_res(i,:),   col.crew,   'Crew weight');
    s_pitch = pitch_sail_res(i,:) + pitch_main_res(i,:) + pitch_rudder_res(i,:) ...
            + pitch_hull_res(i,:) + pitch_crew_res(i,:);
    plot(TWA_vec, s_pitch, '--', 'Color', col.sum, 'LineWidth', 1.4, 'DisplayName', 'Sum (should be ~0)');
    fmt('Pitching moments (My about boat CG, + = bow up)', 'Moment [Nm]');
end


%%
%--------------------------------------------------------------------------
% Crew position plot (top view): wing in grey, crew path coloured by TWA
%--------------------------------------------------------------------------
figure('Name','Crew position','Color','w');
nc = ceil(sqrt(nTWS));  nr = ceil(nTWS/nc);
tl = tiledlayout(nr, nc, 'TileSpacing','compact', 'Padding','compact');

yw   = moth.y_wing_edge;                 % real wing half-span (about 1.09 m)
ylim_crew = -moth.y_crew_foiling;        % crew limit (windward side)
yext = max(yw, ylim_crew) + 0.2;         % axis half-width, leaves room outside the wing
x0 = moth.wing_x_pos;
x1 = x0 + moth.wing_length;

for i = 1:nTWS
    nexttile; hold on; axis equal; box on;

    % Wing area in grey (true wing width)
    patch([-yw yw yw -yw], [x0 x0 x1 x1], [0.80 0.80 0.80], ...
          'EdgeColor',[0.4 0.4 0.4], 'LineWidth',1);
    plot([0 0],[x0 x1],'k:');                                      % centreline

    % Region the solver allows the crew to move in
    rectangle('Position',[-ylim_crew, x0, ylim_crew, moth.wing_length], ...
              'EdgeColor',[0.3 0.3 0.3],'LineStyle','--');

    % Trajectory coloured by TWA
    yy = y_crew_res(i,:);  xx = x_crew_res(i,:);  cc = TWA_vec;
    surface([yy;yy], [xx;xx], zeros(2,numel(yy)), [cc;cc], ...
            'FaceColor','none','EdgeColor','interp','LineWidth',4);

    colormap(turbo); clim([min(TWA_vec) max(TWA_vec)]);
    xlim([-yext yext]);  ylim([x0-0.1 x1+0.1]);
    xlabel('y [m]  (- = windward)'); ylabel('x [m]  (bow up)');
    title(sprintf('TWS = %g m/s', TWS_vec(i)));
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'TWA [deg]';
title(tl, 'Crew position relative to the wing (top view)');


function set_axes(ttl, yl, TWA_vec)
    ax = gca;
    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    grid(ax, 'on');
    set(ax, 'Box', 'on', 'FontSize', 10, 'GridAlpha', 0.15, 'LineWidth', 0.8, ...
            'XLim', [min(TWA_vec) max(TWA_vec)]);
    xlabel(ax, 'TWA [deg]');
    ylabel(ax, yl);
    title(ax, ttl, 'FontWeight', 'bold');
    lg = legend(ax, 'Location', 'eastoutside');
    lg.Box = 'off';
    lg.FontSize = 9;
end

%%
% =========================================================================
% EXTRA PLOTS: Sail Aerodynamics vs AWA (One Subplot per TWS, Single Y-Axis)
% =========================================================================

% Convert all saved AWA from radians to degrees
awa_deg_all = rad2deg(sail_awa_res);

% Find the absolute min and max of the used AWAs to tightly bound the x-axis
min_awa = min(awa_deg_all(:));
max_awa = max(awa_deg_all(:));

% Create a single large figure
fig_aero = figure('Name', 'Sail Aerodynamics Grid', 'Color', 'w', 'Position', [100 100 1200 800]);

% 'flow' automatically calculates the best row/column layout for your number of TWS
tl = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Sail Aerodynamic Coefficients vs Apparent Wind Angle', 'FontWeight', 'bold', 'FontSize', 12);

for k = 1:nTWS
    ax = nexttile; 
    hold(ax, 'on');
    
    awa_deg = awa_deg_all(k,:);
    
    % --- Plot Lift Coefficient (CL) ---
    plot(ax, awa_deg, sail_CL_res(k,:), 'o-', 'Color', cmap(k,:), ...
         'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'C_L (Lift)');
    
    % --- Plot Drag Coefficient (CD) ---
    plot(ax, awa_deg, sail_CD_res(k,:), 's--', 'Color', cmap(k,:), ...
         'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'C_D (Drag)');
    
    % Format subplot
    grid(ax, 'on');
    set(ax, 'Box', 'on', 'FontSize', 9, 'GridAlpha', 0.15, ...
        'XLim', [min_awa max_awa], 'YLim', [0 1.6]);
    
    % Labels and titles
    ylabel(ax, 'Coefficient');
    xlabel(ax, 'AWA [deg]');
    title(ax, sprintf('TWS = %g m/s', TWS_vec(k)), 'Color', cmap(k,:), 'FontWeight', 'bold');
    
    % Add a legend to the first subplot only to avoid cluttering every panel
    if k == 1
        legend(ax, 'Location', 'best');
    end
end