function [y_crew, alpha_sail] = get_controls(s, moth)
% converts the controlvariable s into where the sailor sits sideways
% (y_crew) and how far the sail is eased (alpha_sail).
% Example: s=0.5 puts crew halfway out with the sail at full power
%          s=1.5 has crew fully hiked with the sail halfway between min and
%                max sail angle

s1 = min(max(s,-1),1);   % hiking fraction
s2 = min(max(s-1,0),1); % sail depowering fraction

% linear interpolation between centreline position (0) and the full-hike
% position
y_crew     = moth.y_crew_mid + s1*(moth.y_crew_foiling - moth.y_crew_mid);

% interpolates 
alpha_sail = moth.alpha_sail_max - s2*(moth.alpha_sail_max - moth.alpha_sail_min);


end