% WaveSurface
%
% Created by Anders Rosén, aro@kth.se.
%
% Function creating a wave surface.
%
% INPUT:
% state.H      : wave height    [m]; [1x1] regular wave; [Nx1] irreg; =0=calm water; 
%      .lambda : wave length    [m]; [1x1] regular wave; [Nx1] irreg; =0=calm water; 
%      .eps    : wave phase     [-]; [1x1] regular wave; [Nx1] irreg; =0=calm water; 
%      .my     : wave direction [rad]
%      .t      : time           [s]
% Xsurf        : x coordinates  [m]
% Ysurf        : y coordinates  [m]
%
% OUTPUT:
% Zsurf        : wave level in x,y [m]

function [Zsurf] = WaveSurface(state,Xsurf,Ysurf)

if state.H==0
    Zsurf = zeros(size(Xsurf));
else
    k = 2*pi./(state.lambda + 1e-4);
    w = sqrt(9.81*k);
    Zsurf = zeros(size(Xsurf));
    for I = 1:length(state.H)
        Zsurf = [Zsurf + (state.H(I)/2) * cos(k(I)*(Xsurf*cos(state.my)+Ysurf*sin(state.my)) - w(I)*state.t + state.eps(I))];
    end
end

Zsurf = round(Zsurf*1000)/1000; Zsurf = Zsurf + 1e-4;           % for avoiding bugs

% state.Xsurf = round(state.Xsurf*1000)/1000; state.Ysurf = round(state.Ysurf*1000)/1000; % for avoiding bugs
