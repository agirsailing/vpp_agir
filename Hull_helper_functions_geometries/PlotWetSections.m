% PlotWetSections
%
% Function that plots a hull geometry defined in the structured array "ship"
% which has been aquired by calling the function
%   [ship]=ReadHullGeometry;
% and water lines, wet section areas and the wet section areas' centres of buoyancy
% defined in the structured array "state"
% and that have been determined by calling the function
%   [state] = WetSections(ship,state)

%%
function [] = PlotWetSections(ship,state,WS)

figure(101); pause(0.01); hold off; clf;

% --- plot hull geometry & wet section areas --------------------------
for I=1:ship.ns
    subplot(1,4,1);   % ship fixed coord syst
    plot3(ship.x(1:ship.np(I),I),ship.y(1:ship.np(I),I),ship.z(1:ship.np(I),I),'k-','LineWidth',1); hold on; % geometry
    subplot(1,4,2:4); % global coord syst
    [X,Y,Z] = TransShipfixedGlobal(ship.x,ship.y,ship.z,state.eta);
    plot3(X(1:ship.np(I),I),Y(1:ship.np(I),I),Z(1:ship.np(I),I),'k-','LineWidth',1); hold on; % geometry
    if WS.npw(I)>0
        subplot(1,4,1);   % ship fixed coord syst
        p1 = patch(WS.xw(1:WS.npw(I),I),WS.yw(1:WS.npw(I),I),WS.zw(1:WS.npw(I),I),'c');
        set(p1,'facecolor','cyan','edgecolor','none'); alpha(p1,0.2); %camlight; % wet section areas
        plot3(WS.xwl',WS.ywl',WS.zwl','c.','MarkerSize',8);  %water line
        subplot(1,4,2:4); % global coord syst
        [Xw,Yw,Zw]    = TransShipfixedGlobal(WS.xw,WS.yw,WS.zw,state.eta);
        p2 = patch(Xw(1:WS.npw(I),I),Yw(1:WS.npw(I),I),Zw(1:WS.npw(I),I),'c');
        set(p2,'facecolor','cyan','edgecolor','none'); alpha(p2,0.2); %camlight; % wet section areas
        [Xwl,Ywl,Zwl] = TransShipfixedGlobal(WS.xwl,WS.ywl,WS.zwl,state.eta); %water line
        plot3(Xwl',Ywl',Zwl','c.','MarkerSize',8);
    end
end
subplot(1,4,1); camlight;
subplot(1,4,2:4); camlight;
% ---------------------------------------------------------------------

% --- plot water surface ----------------------------------------------
subplot(1,4,2:4); % global coord syst
maxX=max(max(X)); minX=min(min(X)); maxY=max(max(Y)); minY=min(min(Y));
dX=maxX-minX; intX=(minX-dX*0.2):(dX/100):(maxX+dX*0.2); 
dY=maxY-minY; intY=(minY-dY*0.4):(dY/100):(maxY+dY*0.4);
[Xsurf,Ysurf]=meshgrid(intX,intY);
[Zsurf] = WaveSurface(state,Xsurf,Ysurf);
% Zsurf = zeros(size(Xsurf));
axis equal; axis tight;
p3 = surf(Xsurf,Ysurf,Zsurf); set(p3,'facecolor','cyan','edgecolor','none'); alpha(p3,0.1); %camlight;
% ---------------------------------------------------------------------

% --- plot (xB,yB,zB) -------------------------------------------------
subplot(1,4,1);   % ship fixed coord syst
plot3(WS.xB,WS.yB,WS.zB,'b.');
subplot(1,4,2:4); % global coord syst
[XB,YB,ZB] = TransShipfixedGlobal(WS.xB,WS.yB,WS.zB,state.eta);
plot3(XB,YB,ZB,'b.');
% ---------------------------------------------------------------------

% --- titles, labels, axes --------------------------------------------
subplot(1,4,1); title('Ship fixed coordinate system'); xlabel('x'); ylabel('y'); zlabel('z');
view([-90 0]); axis equal; axis tight; grid on; box on; q1=axis;q1(6)=q1(6)*1.05; axis(q1);
subplot(1,4,2:4);title('Global coordinate system'); xlabel('X'); ylabel('Y'); zlabel('Z'); 
axis equal; axis tight; grid on;
% ---------------------------------------------------------------------
