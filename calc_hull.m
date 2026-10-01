function [R V] = calc_hull(VS,H,hulldata)
% Resistance and displaced volume of the hull, interpolated from the table
% made by make_preCalcHull. H = 0 hull out of the water, H < 0 immersed.

load(hulldata.preCalc); % R_ (rows H_, columns VS_), V_, H_, VS_

if H<0 && H>=H_(1)
    R_vec = interp1(VS_,R_',VS,'linear','extrap');
    R = interp1(H_,R_vec',H,'linear','extrap');
    V = interp1(H_,V_,H,'linear','extrap');
else
    V = 0;
    R = 0;
end

R = max([R 0]);