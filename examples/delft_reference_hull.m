function [hull, env] = delft_reference_hull()
%DELFT_REFERENCE_HULL Published DSYHS 25 geometry at LWL = 10 m.
% de Baar et al. (2015), doi:10.1016/j.compfluid.2015.10.004,
% accepted manuscript Table 2 (p.28): LWL=2, BWL=.5, draft=.0928,
% volume=.03698, Swet=.671, scaled geometrically by 5 below.
% Keuning-Katgert (2008), Table 1, Sysser 25: Cp=.548, Cm=.727,
% LCB_fpp/LWL=.520, LCB_fpp/LCF_fpp=.936, volume^(2/3)/Awp=.165.
% Table ratios are rounded; this is a documented reconstruction, not CAD.
% xFPP=12 m chooses an illustrative stern-deck origin 2 m aft of the aft
% waterline endpoint. Only coordinate differences affect the prediction.
% Water values are explicit example assumptions, not reported tank values.
hull.LWL = 2 * 5;
hull.BWL = 0.5 * 5;
hull.draft = 0.0928 * 5;
hull.volume = 0.03698 * 5^3;
hull.Swet = 0.671 * 5^2;
hull.Awp = hull.volume^(2/3) / 0.165;
hull.Cp = 0.548;
hull.Cm = 0.727;
hull.xFPP = 12;
hull.xLCB = hull.xFPP - 0.520 * hull.LWL;
hull.xLCF = hull.xFPP - 0.520 * hull.LWL / 0.936;
env.g = 9.81;
env.water.rho = 1025;
env.water.nu = 1.19e-6;
end
