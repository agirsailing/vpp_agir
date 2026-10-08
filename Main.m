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
TWS_vec = [5,10];      % [m/s] True wind speeds
TWA_vec = 45:5:155;         % [deg] True wind angles (10-degree steps)

nTWS = length(TWS_vec);
nTWA = length(TWA_vec);

% Pre-allocate results matrices
VS_res = zeros(nTWS, nTWA);
TW_res = zeros(nTWS, nTWA);

% Run the sweeps
for i = 1:nTWS
    TWS = TWS_vec(i);
    fprintf('\n====================================\n');
    fprintf('Starting Wind Speed: %d m/s\n', TWS);
    
    % Initial guess for 20 degrees (using 1.2x TWS because 20 degrees is very tight/slow)
    % 5-DOF initial guess: [VS; heave; heel; pitch; leeway]
    X0 = [TWS * 1.3; moth.depth_foil_flying; deg2rad(0); deg2rad(1.0); deg2rad(3.0)];
        
    for j = 1:nTWA
        TWA_deg = TWA_vec(j);
        TWA_rad = deg2rad(TWA_deg); 
        
        fprintf('  TWA %3d deg... ', TWA_deg);
        
        % Call the 5FOF Newton solver
        [X, iter, FLAG] = solve_Newton(X0,env,moth,TW_sail,TWS,TWA_rad);
        
        VS = X(1);

        % Safety Check: Did the point fail, stall, or hit the 50 iteration limit?
        if FLAG ~= 0 || isnan(VS) || VS < 1 || (moth.mv_span - X(2)) < -0.16
            fprintf('Stalled/Failed. Resetting.\n');
            VS_res(i,j) = NaN; % Record as 0 to indicate it cannot foil here

        else
            fprintf('Success! VS = %5.2f m/s (iters: %d)\n', VS, iter);
            VS_res(i,j) = VS;
            
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

polarplot(TWA_plot_rad, VS_res(1,:), 'b-o', 'LineWidth', 1.5, 'DisplayName', "TWS ="+num2str(TWS_vec(1))+" m/s"); 
hold on;
polarplot(TWA_plot_rad, VS_res(2,:), 'r-o', 'LineWidth', 1.5, 'DisplayName', "TWS ="+num2str(TWS_vec(2))+" m/s");
%polarplot(TWA_plot_rad, VS_res(3,:), 'g-o', 'LineWidth', 1.5, 'DisplayName', 'TWS = 20 m/s');

% Customize the polar axes to match sailing conventions
pax = gca;
pax.ThetaZeroLocation = 'top';      % Places 0 degrees (wind direction) at the top
pax.ThetaDir = 'clockwise';         % Angles increase clockwise
% pax.ThetaAngleFormat = 'degrees';   % Displays the axis labels in degrees

title('Velocity Prediction Polar', 'FontSize', 12);
legend('Location', 'best');