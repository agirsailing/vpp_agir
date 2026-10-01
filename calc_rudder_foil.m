function [diff_residuals,FH, Moments ,CPx,CPy,CPz,alfa_vec] = ...
    calc_foil(surge, heave, heel, pitch, leeway, rake, rudder_angle, sections_yzc,env,moth,plotflag,Target_Loads)
%
% surge        [m/s] Forward boat speed
% heave        [m]   Vertical displacement (positive downwards)
% heel         [rad] Roll angle (positive = starboard down)
% pitch        [rad] Pitch angle (positive = bow up)
% leeway       [rad] Drift angle (positive = drift to starboard)
% rake         [rad] Foil sweep (positive = blade tip aft)
% rudder_angle [rad] Helm angle (positive = trailing edge to starboard)
%
% Optional arguments:
%   plotflag: determines plotting true/false  (1/0)
%   FyFz_target [N] contains vector of forces which is subtracted 
%     from the resulting FH for the objective function for optimizer to solve.
%
% Notation:
%   positive leeway = positive z-rotation
%   positive rake   = positive y-rotation
%   positive heel   = positive x-rotation
%
% Assumptions:
% - Entire foil is submerged
% - All foil panels have zero built-in angle of attack.
%   i.e. all panels are built in the x-z-plane.
% - Simple airfoil in terms on Cl & Cd (Naca 0012)
% - Panel area calcs are coarse ;-)
% - Panel lift is located in panel center.
%
%-------------------------------------------------------------------------
if nargin<11;plotflag = false;end    % To make this argument optional
if nargin<12;Target_Loads = [0;0;0;0;0;0];end % To make this argument optional

thisDir = fileparts(which('calc_foil'));
addpath(genpath(fullfile(thisDir, 'Helper_functions')));
    
% For clarity
q        = 0.5*1000*surge^2;  % [N/m2] Dynamic pressure 

% Cut wet sections, keep only the submerged part of the foil
[Ywet,Zwet] = cut_wet_sections(sections_yzc(1,:),sections_yzc(2,:), heave);
sections_yzc = [Ywet; Zwet; sections_yzc(3,end-length(Ywet)+1:end)];
sections_yzc = sections_yzc';
npanels  = length(sections_yzc(:,1))-1; % Number of wet panels in foil

% Boat Frame Rotation
c1=cos(heel);   s1=sin(heel);   % Pre-calculate the trigonometric functions
c2=cos(pitch);   s2=sin(pitch);   % Pre-calculate the trigonometric functions
c3=cos(leeway); s3=sin(leeway); % Pre-calculate the trigonometric functions

Tx = [  1   0   0;     0 c1 -s1;     0 s1 c1]; % partial Rotation around x-axis
Ty_pitch = [  c2  0  s2;     0  1   0;   -s2  0 c2]; % partial Rotation around y-axis
Tz_leeway = [  c3 -s3  0;    s3 c3   0;     0  0  1]; % partial Rotation around z-axis

% The local angles for the rudders

c4 = cos(rake);         s4 = sin(rake);
c5 = cos(rudder_angle) ; s5 = sin(rudder_angle);

Ty_rake =   [ c4 0 s4 ; 0  1  0; -s4 0 c4];
Tz_rudder = [c5 -s5 0 ; s5 c5 0;  0  0  1];

T  = Tx*Ty_pitch*Tz_leeway*Ty_rake*Tz_rudder; % Total transformation matrix (order matters!!)

% Now we need to do some panel-geometry calculations.  Loop through your 
% panels and rotate them according to Heel, leeway and rake.
% Notation: North is up, west is forward. 
foil_length = 0;   % [m] Just  initialization
for ipanel=1:npanels
  % Calc the panel corner coords
  %                      x                            y                      z
  nw     = T*[ sections_yzc(ipanel,3)/2  ; sections_yzc(ipanel,1)  ;sections_yzc(ipanel,2)  ];
  ne     = T*[-sections_yzc(ipanel,3)/2  ; sections_yzc(ipanel,1)  ;sections_yzc(ipanel,2)  ];
  sw     = T*[ sections_yzc(ipanel+1,3)/2; sections_yzc(ipanel+1,1);sections_yzc(ipanel+1,2)];
  se     = T*[-sections_yzc(ipanel+1,3)/2; sections_yzc(ipanel+1,1);sections_yzc(ipanel+1,2)];
  mid    = (nw+ne+sw+se)/4;                          % [m]   Panel mid point
  ni     = cross((nw-ne),(sw-nw));ni=ni/norm(ni);    % [-]   Panel normal
  alfai  = (-pi/2+acos(dot([1,0,0],ni)));            % [rad] Panel angle of attack
  Ai     = norm(cross((nw-ne),(sw-nw)));             % [m2]  Panel area
  % Collect results into vectors
  A_vec(ipanel)    = Ai;        % Collect all panel areas in a vector
  alfa_vec(ipanel) = alfai;     % Collect all alfas in a vector
  Ni_vec(:,ipanel) = ni;        % Collect all normals in a vector
  mid_vec(:,ipanel)= mid;       % Collect all pidpoints in a "vector"
  foil_length      = foil_length+norm(nw(2:3)-se(2:3)); % sum up the foil length
  % No do some Plotting
  if plotflag;
    %plot3([nw(1) ne(1) se(1) sw(1) nw(1)],[nw(2) ne(2) se(2) sw(2) nw(2)],[nw(3) ne(3) se(3) sw(3) nw(3)]);hold on;
    patch([nw(1) ne(1) se(1) sw(1)],[nw(2) ne(2) se(2) sw(2)],[nw(3) ne(3) se(3) sw(3)],[1 1 1 1],'FaceAlpha',0.4);hold on;
    plot3([mid(1)-0 mid(1)+1],[mid(2) mid(2)],[mid(3) mid(3)],'g'); % Free flow vector
    plot3([mid(1) mid(1)+ni(1)],[mid(2) mid(2)+ni(2)],[mid(3) mid(3)+ni(3)],'g');
  end
end
AR = foil_length^2/sum(A_vec); % [-] Aspect ratio
e  = 0.85;                     % [-] Span efficiency factor
chord = moth.f2_c_root;
nu = env.water.nu;
Re = surge*moth.f2_c_root/env.water.nu; % reynolds number
% Now, most of the geometric pre-calculations are done... puhh

foil_data = moth.r_foildata;
% extract data corresponding to the closest Re from table
Re_table = extract_data(foil_data,Re);

% make alpha cl and alpha cd array 
a_cl = [Re_table.alpha,Re_table.CL];
a_cd = [Re_table.alpha,Re_table.CD];

FH=[0,0,0]';
% test for interpolation in data table instead of empirical formulas
for ipanel=1:npanels
    alfa_deg = alfa_vec(ipanel) * (180/pi);
    Cl   = interp1(a_cl(:,1),a_cl(:,2),alfa_deg);    % [-] 2D-lift interpolated from xfoil data.
    Cd   = interp1(a_cd(:,1),a_cd(:,2),alfa_deg);    % [-] 2D-drag interpolated from xfoil data.
    CL   = Cl/(1+2/(e*AR));             % [-] 3D lift-correction.
    L    = q*CL*A_vec(ipanel);          % [N] Total panel lift (in the yz-plane).
    CDi  = CL^2/(pi*e*AR);              % [-] Induced drag coeff.
    D    = q*A_vec(ipanel)*(Cd+CDi);    % [N] Total drag.
    Ni   = Ni_vec(:,ipanel);            % [-] Panel normal.
    niyz = [0;Ni(2);Ni(3)];niyz = niyz/norm(niyz); % Unit panel normal projection on the yz-plane.
    FHi  = niyz*L +[-D;0;0];            % Panel force vector.
    FH_panels(:,ipanel)  = FHi;         % Array of panel forces.
end


% for ipanel=1:npanels
%   Cl_a = 2*pi;                        % [-] 2D Lift derivative for ?NACA0012.   
%   Cl   = Cl_a*alfa_vec(ipanel);       % [-] 2D lift.
%   Cd   = 0.008 + 0.01*Cl^2;           % [-] 2D-drag for ?NACA0012.
%   CL   = Cl/(1+2/(e*AR));             % [-] 3D lift-correction.
%   L    = q*CL*A_vec(ipanel);          % [N] Total panel lift (in the yz-plane).
%   CDi  = CL^2/(pi*e*AR);              % [-] Induced drag coeff.
%   D    = q*A_vec(ipanel)*(Cd+CDi);    % [N] Total drag.
%   Ni   = Ni_vec(:,ipanel);            % [-] Panel normal.
%   niyz = [0;Ni(2);Ni(3)];niyz = niyz/norm(niyz); % Unit panel normal projection on the yz-plane. 
%   FHi  = niyz*L +[-D;0;0];            % Panel force vector.
%   FH_panels(:,ipanel)  = FHi;         % Array of panel forces.
% end


% Now sum upp each panels contribution into FH.
if npanels>1;FH = sum(FH_panels')'; else FH=FH_panels;end      % Total hydrodynamic force from all panels



% Calc the moments around attachment point
Fx_ = FH_panels(1,:);  % [N] Panels x-contributions
Fy_ = FH_panels(2,:);  % [N] All panel contributions-vector
Fz_ = FH_panels(3,:);  % [N] All panel contributions-vector
x_  = mid_vec(1,:);    % [m] Panels mid-points
y_  = mid_vec(2,:);    % [m] Panels mid-points
z_  = mid_vec(3,:);    % [m] Panels mid-points
Moments = [Fz_*y_'+Fy_*z_';Fx_*z_'+Fz_*x_';Fx_*y_'+Fy_*x_']; % [Nm] Total moments around the 3 axes at clamping
 
% Target_Loads = [Fx_target; Fy_target; Fz_target; Mx_target; My_target; Mz_target]
actual_loads = [FH(1); FH(2); FH(3); Moments(1); Moments(2); Moments(3)];
diff_residuals = Target_Loads - actual_loads;

% Potential problem we need to deal with...
% When forces are zero => CP is undefined :-(
% Ugly fix:
if abs(FH(1))  <0.00000001;FH(1)  =0.00000001;end  % Needed for CP-calcs
if abs(FH(2))  <0.00000001;FH(2)  =0.00000001;end  % Needed for CP-calcs
if abs(FH(3))  <0.00000001;FH(3)  =0.00000001;end  % Needed for CP-calcs

% Now calculate center of pressure for entire foil.
CPx = [0;sum(Fx_.*y_);sum(Fx_.*z_)]/FH(1);% [m] Centre of pressure in x-dir
CPy = [sum(Fy_.*x_);0;sum(Fy_.*z_)]/FH(2);% [m] Centre of pressure in y-dir
CPz = [sum(Fz_.*x_);sum(Fz_.*y_);0]/FH(3);% [m] Centre of pressure in z-dir

if plotflag;
 for ipanel=1:npanels;
   nFHi = FH_panels(:,ipanel)*2/max(norm(FH_panels(:,:))); % Init vector in the der of the panel force
   midi = mid_vec(:,ipanel);
   plot3([midi(1),midi(1)+nFHi(1)],[midi(2),midi(2)+nFHi(2)],[midi(3),midi(3)+nFHi(3)],'r','Linewidth',2);
 end
 % Plot the lines for centre of pressure
 plot3([-100;100],[CPx(2),CPx(2)],[CPx(3) CPx(3)],'Color',[0.5 0.5 0.5]);
 plot3([CPy(1),CPy(1)],[-100;100],[CPy(3) CPy(3)],'Color',[0.5 0.5 0.5]);
 plot3([CPz(1),CPz(1)],[CPz(2) CPz(2)],[-100;100],'Color',[0.5 0.5 0.5]);
 rotate3d on;grid on;
 xlabel('X');ylabel('Y');zlabel('Z')
 set(gca, 'ZDir', 'reverse'); % Visually flips Z-axis to point downwards
 axis([-2 4 -1 5 0 6]);       % Z-limits updated for positive depth

 fprintf('---------------------------------------------------\n');
 fprintf('FH = [%.0f,%.0f,%.0f] N   => Total=%.0fN \n',FH(1),FH(2),FH(3),norm(FH));
 fprintf('Foil area           %.2f m2  \n',sum(A_vec));
 fprintf('Foil aspect ratio   %.1f m2  \n',AR);
 fprintf('Max angle of attack %.1f deg \n',max(alfa_vec)*180/pi);
 fprintf('Min angle of attack %.1f deg \n',min(alfa_vec)*180/pi);
end


fprintf('Surge  = %.2f m/s\n', surge);
fprintf('Heave  = %.2f m\n', heave);
fprintf('Heel   = %.1f deg\n', rad2deg(heel));
fprintf('Pitch  = %.1f deg\n', rad2deg(pitch));
fprintf('Leeway = %.1f deg\n', rad2deg(leeway));
fprintf('Rake   = %.1f deg\n', rad2deg(rake));
fprintf('Rudder = %.1f deg\n', rad2deg(rudder_angle));