function [y_crew, alpha_sail] = get_controls(s, moth)
% s in [0,2]: 0..1 crew hikes out, 1..2 sail depowers
s1 = min(max(s,0),1);
s2 = min(max(s-1,0),1);
y_crew     = moth.y_crew_mid + s1*(moth.y_crew_foiling - moth.y_crew_mid);
alpha_sail = moth.alpha_sail_max - s2*(moth.alpha_sail_max - moth.alpha_sail_min);
end