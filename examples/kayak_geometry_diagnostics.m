function audit = kayak_geometry_diagnostics(result)
%KAYAK_GEOMETRY_DIAGNOSTICS Inspect an existing kayak_resistance result.
% Sampled section maximum is diagnostic only; resistance coefficients unchanged.
hs=result.hydrostatics; mesh=result.mesh;
x=linspace(hs.waterlineExtent(1),hs.waterlineExtent(2),201);
% End caps have nonzero area: sample just inside both endpoints.
dx=1e-7*hs.LWL; x(1)=x(1)+dx; x(end)=x(end)-dx;
area=arrayfun(@(s)hgeom.immersed_section_area(hs.mesh,s),x);
[maximum,k]=max(area);
audit=struct('x',x,'area',area,'sampledMaximumArea',maximum, ...
    'sampledMaximumX',x(k),'closureAdjusted',sum([mesh.sections.closureAdjusted]));
fig=figure('Name','Provisional kayak geometry','Position',[100 100 1100 750]);
audit.figure=fig; layout=tiledlayout(fig,2,2,'Padding','loose','TileSpacing','loose');
v=mesh.vertices-[0 0 result.floating.pose.elevation];
ax=nexttile; trisurf(mesh.faces,v(:,1),v(:,2),v(:,3),'EdgeColor','none','FaceAlpha',.8);
axis(ax,'equal'); set(ax,'ZDir','reverse'); view(ax,35,25);
xlabel(ax,'x forward [m]'); ylabel(ax,'y starboard [m]'); zlabel(ax,'Depth [m]');
title(ax,sprintf('Sealed mesh: %d adjusted sections',audit.closureAdjusted));
ax=nexttile; trisurf(mesh.faces,v(:,1),v(:,2),v(:,3),'EdgeColor','none');
view(ax,0,0); axis(ax,'equal'); set(ax,'ZDir','reverse'); hold(ax,'on');
plot3(ax,hs.waterlineExtent,repmat(min(v(:,2))-.01,1,2),[0 0],'r-','LineWidth',2);
xlabel(ax,'x [m]'); zlabel(ax,'Depth [m]'); title(ax,'Profile and waterline (red)');
ax=nexttile; trisurf(hs.mesh.faces,hs.mesh.vertices(:,1),hs.mesh.vertices(:,2),hs.mesh.vertices(:,3),'EdgeColor','none');
view(ax,2); axis(ax,'equal'); xlabel(ax,'x [m]'); ylabel(ax,'y [m]'); title(ax,'Immersed mesh, plan view');
ax=nexttile; plot(ax,x,area,'LineWidth',1.5); hold(ax,'on');
plot(ax,result.hull.midshipX,result.hull.midshipArea,'o',x(k),maximum,'x');
xlabel(ax,'x [m]'); ylabel(ax,'Immersed section area [m^2]'); grid(ax,'on');
legend(ax,{'Sampled sections','Current midship','Sampled maximum'},'Location','best');
layout.Units='normalized'; layout.OuterPosition=[0 0 1 .91];
annotation(fig,'textbox',[0 .94 1 .05], ...
    'String','Provisional kayak geometry — unverified import and sealed closures', ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',13);
set(findall(fig,'Type','axes'),'FontSize',9);
end
