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
TW_sail = deg2rad(10); % [rad] Assumed sail twist

takeoff_speed = 5; %[m/s]
%--------------------------------------------------------------------------    
% Enviornmental variables, Loading
%--------------------------------------------------------------------------    

%--------------------------------------------------------------------------
% BEGIN VPP
%--------------------------------------------------------------------------


% % Exercise 1 here...
% calc_residuals_Newton(VS,TW,H,TWS,TWA,hulldata,wingdata);
% VS0=14.4;
% TW0=0.1745;
% [VS,TW,iter,FLAG] = solve_Newton(VS0,TW0,H,TWS,TWA,hulldata,wingdata);


%--------------------------------------------------------------------------
% BEGIN VPP - Final Polar Sweep (20 to 160 degrees)
%--------------------------------------------------------------------------
H = -0.5; % [m] Constant flying height

% Define vectors for the required homework span
TWS_vec = [1,2,3,4,5,6,7,8];      % [m/s] True wind speeds
TWA_vec = 40:5:155;         % [deg] True wind angles (10-degree steps)

nTWS = length(TWS_vec);
nTWA = length(TWA_vec);

% Pre-allocate results matrices
VS_res     = NaN(nTWS, nTWA);   % surge speed [m/s]
heave_res  = NaN(nTWS, nTWA);   % heave [m]
heel_res   = NaN(nTWS, nTWA);   % heel [deg]
pitch_res  = NaN(nTWS, nTWA);   % pitch [deg]
leeway_res = NaN(nTWS, nTWA);   % leeway [deg]

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

% Pitch moments [Nm]
pitch_sail_res    = NaN(nTWS, nTWA);
pitch_hull_res    = NaN(nTWS, nTWA);
pitch_main_res    = NaN(nTWS, nTWA);
pitch_rudder_res  = NaN(nTWS, nTWA);

% Run the sweeps
for i = 1:nTWS
    TWS = TWS_vec(i);
    fprintf('\n====================================\n');
    fprintf('Starting Wind Speed: %d m/s\n', TWS);

    if TWS < takeoff_speed
        % Displacement mode: hull immersed (~10-15 cm draught)
        heave_0 = moth.mv_span + 0.12; % ~1.12 m
        VS_0    = TWS * 0.8;
    else
        % Foiling mode: main foil flying near target depth
        heave_0 = moth.depth_foil_flying; % 0.18 m
        VS_0    = TWS * 1.5;
    end

    X0 = [VS_0; heave_0; deg2rad(0); deg2rad(1.0); deg2rad(-3.0)];

    % % 5-DOF initial guess: [VS; heave; heel; pitch; leeway]
    % X0 = [TWS * 1.3; 0.88; deg2rad(0); deg2rad(1.0); deg2rad(-3.0)];
    % 
    % 
    for j = 1:nTWA
        TWA_deg = TWA_vec(j);
        TWA_rad = deg2rad(TWA_deg); 
        
        fprintf('  TWA %3d deg... ', TWA_deg);
        
        % Call the 5FOF Newton solver
        [X, iter, FLAG, FM] = solve_Newton(X0,env,moth,TW_sail,TWS,TWA_rad);
        disp(" ")
        disp(FM)
        VS = X(1);
        heave = X(2);
        heel = rad2deg(X(3));
        pitch = rad2deg(X(4));
        leeway = rad2deg(X(5));

        % Safety Check: Did the point fail, stall, or hit the 50 iteration limit?
        if FLAG ~= 0 || isnan(VS) || VS < 1 || (moth.mv_span - X(2)) < -0.16
            fprintf('Stalled/Failed. Resetting.\n');
            disp("Surge="+num2str(VS)+"m/s, Heave="+num2str(heave)+"m, Heel="+num2str(heel)+"deg, Pitch="+num2str(pitch)+"deg, leeway="+num2str(leeway)+"deg.")
            VS_res(i,j) = NaN; % Record as 0 to indicate it cannot foil here
            heave_res(i,j)  = NaN;
            heel_res(i,j)   = NaN;
            pitch_res(i,j)  = NaN;
            leeway_res(i,j) = NaN;


        else
            fprintf('Success! VS = %5.2f m/s (iters: %d)\n', VS, iter);
            VS_res(i,j) = VS;

            heave_res(i,j)  = heave;
            heel_res(i,j)   = heel;
            pitch_res(i,j)  = pitch;
            leeway_res(i,j) = leeway;

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

            % Pitch components
            pitch_sail_res(i,j)   = FM.pitch_sail;
            pitch_hull_res(i,j)   = FM.pitch_hull;
            pitch_main_res(i,j)   = FM.pitch_main_f;
            pitch_rudder_res(i,j) = FM.pitch_rudder_f;

            if X(5) > 0
                fprintf('  (warning: positive leeway, foils drifting to windward)\n');
            end

            disp("Surge="+num2str(VS)+"m/s, Heave="+num2str(heave)+"m, Heel="+num2str(heel)+"deg, Pitch="+num2str(pitch)+"deg, leeway="+num2str(leeway)+"deg.")
            % Update guess for the next angle (Sequential Updating)
            X0 = X;
        end
    end
end


%--------------------------------------------------------------------------
% Plot the Final Results as a Polar Plot
%--------------------------------------------------------------------------
%%
figure;

% polarplot requires the angle array to be in radians
TWA_plot_rad = deg2rad(TWA_vec);

k = 0;
for speed = TWS_vec
    k = k + 1;

    polarplot(TWA_plot_rad, VS_res(k,:), '-o', 'LineWidth', 1.5, 'DisplayName', "TWS ="+num2str(speed)+" m/s");
    hold on
end

grid on
grid minor
% Customize the polar axes to match sailing conventions
pax = gca;
pax.ThetaZeroLocation = 'top';      % Places 0 degrees (wind direction) at the top
pax.ThetaDir = 'clockwise';         % Angles increase clockwise
% pax.ThetaAngleFormat = 'degrees';   % Displays the axis labels in degrees
thetalim([0,180])

title("Velocity Prediction Polar, Take off speed "+num2str(takeoff_speed)+" m/s", 'FontSize', 12);
legend('Location', 'best');
%%
figure;
subplot(5,1,1)
hold on
grid on
grid minor
plot(TWA_vec,VS_res)
ylabel("speed [m/s]")
legend(string(TWS_vec))

subplot(5,1,2)
hold on
grid on
grid minor
plot(TWA_vec,heave_res)
ylabel("heave [m]")

subplot(5,1,3)
hold on
grid on
grid minor
plot(TWA_vec,heel_res)
ylabel("heel [deg]")

subplot(5,1,4)
hold on
grid on
grid minor
plot(TWA_vec,pitch_res)
ylabel("pitch [deg]")

subplot(5,1,5)
hold on
grid on
grid minor
plot(TWA_vec,leeway_res)
ylabel("leeway [deg]")
xlabel("TWA [deg]")

%--------------------------------------------------------------------------
% Plot Force and Moment Balances (Component Breakdown)
%--------------------------------------------------------------------------
for i = 1:nTWS
    curr_tws_kn = TWS_vec(i) * env.ms2kn;
    
    figure('Name', sprintf('Force & Moment Breakdown - TWS %.1f kn', curr_tws_kn), ...
           'Color', 'w', 'Position', [100, 100, 1100, 800]);

    % 1. SURGE BALANCE
    subplot(2,2,1);
    plot(TWA_vec, surge_sail_res(i,:),   'b-',  'LineWidth', 1.8, 'DisplayName', 'Sail Drive (+)');
    hold on;
    plot(TWA_vec, surge_wings_res(i,:),  'm--', 'LineWidth', 1.4, 'DisplayName', 'Wings Drag (-)');
    plot(TWA_vec, surge_main_res(i,:),   'r--', 'LineWidth', 1.6, 'DisplayName', 'Main Foil Drag (-)');
    plot(TWA_vec, surge_rudder_res(i,:), 'k--', 'LineWidth', 1.4, 'DisplayName', 'Rudder Drag (-)');
    plot(TWA_vec, surge_hull_res(i,:),   'c:',  'LineWidth', 1.4, 'DisplayName', 'Hull Drag (-)');
    yline(0, 'k:');
    grid on; grid minor;
    xlabel('TWA [deg]'); ylabel('Force [N]');
    title(sprintf('Surge Forces (X) - TWS %.0f kn', curr_tws_kn));
    legend('Location', 'best');

    % 2. HEAVE BALANCE
    subplot(2,2,2);
    plot(TWA_vec, heave_sail_res(i,:),   'b-',  'LineWidth', 1.6, 'DisplayName', 'Sail Vertical Force');
    hold on;
    plot(TWA_vec, heave_main_res(i,:),   'r-',  'LineWidth', 1.8, 'DisplayName', 'Main Foil Lift');
    plot(TWA_vec, heave_rudder_res(i,:), 'k--', 'LineWidth', 1.4, 'DisplayName', 'Rudder Elevator Lift');
    plot(TWA_vec, heave_hull_res(i,:),   'c:',  'LineWidth', 1.4, 'DisplayName', 'Hull Buoyancy');
    plot(TWA_vec, heave_weight_res(i,:), 'm-',  'LineWidth', 1.4, 'DisplayName', 'Total Weight (+mg)');
    yline(0, 'k:');
    grid on; grid minor;
    xlabel('TWA [deg]'); ylabel('Force [N] (+ = down)');
    title('Heave Forces (Z)');
    legend('Location', 'best');

    % 3. HEEL MOMENT BALANCE
    subplot(2,2,3);
    plot(TWA_vec, heel_sail_res(i,:),   'b-', 'LineWidth', 1.8, 'DisplayName', 'Sail Heel Moment (+)');
    hold on;
    plot(TWA_vec, heel_main_res(i,:),   'r-', 'LineWidth', 1.8, 'DisplayName', 'Main Foil Roll Moment');
    plot(TWA_vec, heel_rudder_res(i,:), 'k-', 'LineWidth', 1.4, 'DisplayName', 'Rudder Roll Moment');
    yline(0, 'k:');
    grid on; grid minor;
    xlabel('TWA [deg]'); ylabel('Moment [Nm]');
    title('Heeling Moments (Mx about CG)');
    legend('Location', 'best');

    % 4. PITCH MOMENT BALANCE
    subplot(2,2,4);
    plot(TWA_vec, pitch_sail_res(i,:),   'b-',  'LineWidth', 1.6, 'DisplayName', 'Sail Pitch Moment');
    hold on;
    plot(TWA_vec, pitch_main_res(i,:),   'r-',  'LineWidth', 1.8, 'DisplayName', 'Main Foil Pitch Moment');
    plot(TWA_vec, pitch_rudder_res(i,:), 'k--', 'LineWidth', 1.6, 'DisplayName', 'Rudder Pitch Moment');
    plot(TWA_vec, pitch_hull_res(i,:),   'c:',  'LineWidth', 1.4, 'DisplayName', 'Hull Pitch Moment');
    yline(0, 'k:');
    grid on; grid minor;
    xlabel('TWA [deg]'); ylabel('Moment [Nm] (+ = bow up)');
    title('Pitching Moments (My about CG)');
    legend('Location', 'best');
end