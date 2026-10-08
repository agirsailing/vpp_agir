function genJava(Re, Foil ,filepath)
%==========================================================================
% Generates input file for JavaFoil
%==========================================================================

% Set options
Country          = '0';                                % Units and such
TransitionModel  = '1';                 % Extended Eppler see reference manual
StallModel       = '0';
                  
% Set analysis conditions
ReFirst    = num2str(round(Re));
ReLast     = num2str(round(Re));
ReDelta    = num2str(round(Re));                % Reynolds number
AngleFirst = '0';          % First angle of attack
AngleLast  = '0';        % Final angle of attack
AngleDelta = '0';       %num2str(Alpha(2)-Alpha(1)); % Step in alpha
TransUpper = '100';               % Default
TransLower = '100';               % Default
Roughness  = '1';                 % Surface roughness, default=smooth=0
AddToPlotsFlag = '0';             % Default

% Create file using filepath
fid = fopen(filepath,'w');

% Header general information
line{1} = (['// recorded on ',datestr(now),' by ',getenv('USER')]);
line{2} = ['Options.Country(',Country,')']; % Settings
line{3} = ['Options.TransitionModel(',TransitionModel,')']; % Settings
line{4} = ['Options.StallModel(',StallModel,')']; % Settings
line{5} = ['Geometry.Open("',pwd,'/JavaFoil/',Foil,'")']; % Generate foil geometry
line{6} = ['Polar.Analyze(',ReFirst,':',ReLast,':',ReDelta,':',AngleFirst,':',AngleLast,':',AngleDelta,':',TransUpper,':',TransLower,':',Roughness,':',AddToPlotsFlag,')']; % Analyze foil
line{7} = ['Polar.Save("',filepath,'")'];
line{8} = ['Exit()'];

% Write file
for i=1:length(line)
    fprintf(fid,'%s\n',line{i});
end

% Close file return to main
fclose(fid);

