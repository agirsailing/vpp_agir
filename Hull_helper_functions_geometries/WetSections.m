% WetSections
%
% Created by Anders Rosén, aro@kth.se. By using this script you agree on:
% - respecting copyright, not spreading the code to anyone else;
% - taking full personal responsibility for how you use the code and
%   results generated.
%
% [WS] = WetSections(ship,state)
%
% Function that determines the intersection between
% a ship hull defined in a ship fixed coordinate system
% and a water surface defined in the global coordinate system
% for a certain hull position and aspect in the global system,
% and the corresponding wet section areas, waterlines, etc
% which can be used for hydrostatic calculations.
%
% INPUT:
% Ship hull geometry:
% ship.x .y .z : hull offset coordinates in a ship fixed coordinate system
%     .ns      : number of sections in the hull geometry definition
%     .np      : number of offset points on the respective sections
%
% Ship position and aspect in the global coordinate system:
% state.eta    : [1x6] vector
%                state.eta(1:3)=[X Y Z] in metres
%                state.eta(4:6)~[roll pitch yaw] in radians
%
% Water surface:
% For calculations in calm water, leave the following variables undefined.
% For calculations in waves:
% state.H      : wave height; [1x1] regular wave; [Nx1] irregular wave; 
%      .lambda : wave length; [1x1] regular wave; [Nx1] irregular wave; 
%      .eps    : wave phase;  [1x1] regular wave; [Nx1] irregular wave;
%      .my     : wave direction
%      .t      : time step         
%
% state.plotflag : =1 plot in local&global coord syst, =2 plot only in global; =0 no plot
%
% OUTPUT:
% WS.Aws                 : wet section areas [m^2]
% WS.xB; WS.yB; WS.zB    : wet section areas' centres of buoyancy in local coordinate system [m]
% WS.Ss                  : wet section arc lengths [m]
% WS.xw; WS.yw; WS.zw    : wet offset points including new waterline in local coordinate system [m]
% WS.npw                 : number of wet offset points including new waterline [-]
% WS.xwl; WS.ywl; WS.zwl : waterline coordinates in local coordinate system [m]
% WS.bwl                 : waterline width at each section [m]
%
% SUB-SCRIPTS:
% TransShipfixedGlobal, TransGlobalShipfixed, TransMat : Euler transformations
% polygonArea, polygonCentroid : for calculating area and area centre of polygons
% WaveSurface : for defining the water surface (calm water or regular waves)

%%
function [WS] = WetSections(ship,state)

% --- initiating variables --------------------------------------------
if isfield(state,'H')        ==0; state.H        =0; end
if isfield(state,'lambda')   ==0; state.lambda   =0; end
if isfield(state,'eps')      ==0; state.eps      =0; end
if isfield(state,'my')       ==0; state.my       =0; end
if isfield(state,'t')        ==0; state.t        =0; end
if isfield(state,'plotflag') ==0; state.plotflag =0; end
if state.H==0; state.lambda=0; state.eps=0;   end
x = ship.x; y = ship.y; z = ship.z; ns = ship.ns; np = ship.np;
Xw   = zeros(size(x,1)*2,ns); Yw = zeros(size(x,1)*2,ns); Zw = zeros(size(x,1)*2,ns);
xw   = []; yw = []; zw  = [];
isw  = zeros(1,ns); npw = zeros(1,ns);
xwl  = zeros(2,ns); ywl = zeros(2,ns); zwl = zeros(2,ns);
iswl = zeros(1,ns); bwl = zeros(1,ns);
Aws  = zeros(1,ns); Ss  = zeros(1,ns);
xB   = zeros(1,ns); yB  = zeros(1,ns); zB = zeros(1,ns); 
% ---------------------------------------------------------------------

% --- transform hull geom fr ship fixed to global coord syst ----------
[X,Y,Z] = TransShipfixedGlobal(x,y,z,state.eta);
% ---------------------------------------------------------------------

% --- determine water surface in X&Y-locations of hull offsets  -------
[Zsurf] = WaveSurface(state,X,Y);
% ---------------------------------------------------------------------

for J = 1:ns % --- stepping through all sections ----------------------
    Tlok  = Zsurf(1:np(J),J) - Z(1:np(J),J); % OP draught, Tlok>0 means wet
    iw    = find(Tlok>0);                    % index of wet OP
    nw    = length(iw);                      % nof wet OP
    
    if nw==0                 % --- whole section dry ------------------
        npw(J) = 0; bwl(J) = 0; Aws(J) = 0; yB(J) = 0; zB(J) = 0;
        
    else
        
        if nw==np(J)         % --- whole section wet ------------------
            isw(J) = J; npw(J) = np(J); bwl(J) = 0;
            Xw(1:npw(J),J)  = X(1:np(J),J); Yw(1:npw(J),J) = Y(1:np(J),J); Zw(1:npw(J),J)=Z(1:np(J),J);
        end
        
        if nw>0 && nw<np(J)  % --- section partly wet -----------------
            isw(J) = J; iswl(J) = J;
            PiPip1x = diff([X(1:np(J),J);X(1,J)]); % x-dist btw offsets Pi and Pi+1
            PiPip1y = diff([Y(1:np(J),J);Y(1,J)]); % y-dist btw offsets Pi and Pi+1
            PiPip1z = diff([Z(1:np(J),J);Z(1,J)]); % z-dist btw offsets Pi and Pi+1
            iwl0 = find(Tlok(2:np(J))./Tlok(1:(np(J)-1))<0); % index offset after WL intersect
            nwl0 = length(iwl0);                   % nof WL intersections
            XT = X(iwl0,J)+PiPip1x(iwl0)./(1+abs(Tlok(iwl0+1)./Tlok(iwl0))); % WL coord
            YT = Y(iwl0,J)+PiPip1y(iwl0)./(1+abs(Tlok(iwl0+1)./Tlok(iwl0))); % WL coord
            ZT = Z(iwl0,J)+PiPip1z(iwl0)./(1+abs(Tlok(iwl0+1)./Tlok(iwl0))); % WL coord
            np0 = np(J);
            Xw0 = X(1:np(J),J); Yw0 = Y(1:np(J),J); Zw0 = Z(1:np(J),J);
            for I=1:nwl0  % adding water line intersection points to original offsets
                Xw0 = [Xw0(1:iwl0(I)) ; XT(I) ; Xw0((iwl0(I)+1):np0)];
                Yw0 = [Yw0(1:iwl0(I)) ; YT(I) ; Yw0((iwl0(I)+1):np0)];
                Zw0 = [Zw0(1:iwl0(I)) ; ZT(I) ; Zw0((iwl0(I)+1):np0)];
                iw  = [iw(iw<=iwl0(I)) ; iwl0(I)+1 ; iw(iw>iwl0(I))+1];
                np0 = np0+1; iwl0 = iwl0+1;
            end
            npw(J) = length(iw);
            Xw(1:npw(J),J) = Xw0(iw); Yw(1:npw(J),J) = Yw0(iw); Zw(1:npw(J),J) = Zw0(iw);
            [xT,yT,zT] = TransGlobalShipfixed(XT,YT,ZT,state.eta);
            [slask,imin] = min(yT); [slask,imax] = max(yT);
            xwl(:,J) = [xT(imin);xT(imax)]; ywl(:,J) = [yT(imin);yT(imax)]; zwl(:,J)=[zT(imin);zT(imax)];
            bwl(J) = sqrt(diff(xwl(:,J)).^2 + diff(ywl(:,J)).^2 + diff(zwl(:,J)).^2);
        end
        
        [xw,yw,zw] = TransGlobalShipfixed(Xw,Yw,Zw,state.eta);
        Aws(J)     = polygonArea(yw(1:npw(J),J),zw(1:npw(J),J));
        yzB        = polygonCentroid(yw(1:npw(J),J),zw(1:npw(J),J));
        yB(J)=yzB(1,1); zB(J)=yzB(1,2); xB=x(1,:);
        if isnan(Aws(J)); Aws(J)=0; end; if isnan(yB(J)); yB(J)=0; end; if isnan(zB(J)); zB(J)=0; end;
        if Aws(J)==0; yB(J)=0; zB(J)=0; end
        sPiPip1    = sqrt(diff(yw(1:npw(J),J)).^2 + diff(zw(1:npw(J),J)).^2); %dist btw wet OP
        Ss(J)      = sum(sPiPip1) - bwl(J); %wet sect arc length
        
    end
    
end % stepping through all sections -----------------------------------

WS.Aws = Aws; WS.xB=xB; WS.yB=yB; WS.zB=zB; WS.bwl=bwl; WS.Ss=Ss;
WS.xw  = xw; WS.yw = yw; WS.zw = zw; WS.npw = npw;
WS.xwl = xwl(:,nonzeros(iswl)); WS.ywl = ywl(:,nonzeros(iswl)); WS.zwl = zwl(:,nonzeros(iswl)); 

if state.plotflag == 1; PlotWetSections(ship,state,WS);              end;
if state.plotflag == 2; PlotWetSectionsOnlyGlobCoord(ship,state,WS); end;
% plotflag=state.plotflag; if exist('plotflag','var'); if state.plotflag == 1; PlotWetSections(ship,state,WS); end; end
% plotflag=state.plotflag; if exist('plotflag','var'); if state.plotflag == 2; PlotWetSectionsOnlyGlobCoord(ship,state,WS); end; end

