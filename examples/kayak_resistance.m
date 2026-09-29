function result = kayak_resistance(totalMass, options)
%KAYAK_RESISTANCE Provisional BRI -> loaded hydrostatics -> Delft check.
% totalMass includes kayak, paddler and gear [kg], default 130.
% Import conventions are UNVERIFIED; all 16 sections need sealed closure.
if nargin<1, totalMass=130; end
if nargin<2, options=struct(); end
root=fileparts(fileparts(mfilename('fullpath')));
importOptions=struct('lengthScale',1,'xHeaderColumn',2,'axisSigns',[1 1 -1], ...
    'originOffset',[0 0 0],'halfHull',true,'halfHullClosure','seal');
mesh=read_bri(fullfile(root,'data','kayak.bri'),importOptions);
env=get_env_params();
% CG is a placeholder required by the solver. Fixed-attitude draft depends
% only on mass; this does NOT solve pitch/roll balance or claim a measured CG.
loading=struct('mass',totalMass,'cg',mean(mesh.vertices,1));
floating=solve_float(mesh,loading,struct('mode','draft','heel',0,'trim',0),env);
assert(floating.converged,floating.message);
[h,hs]=resistance_hull(mesh,floating.pose,env);
V=[0 linspace(.15,.75,121)]*sqrt(env.g*h.LWL);
result=struct('hull',h,'hydrostatics',hs,'floating',floating,'speed',V, ...
    'resistance',[],'details',[],'delftError','', ...
    'importOptions',importOptions,'loading',loading,'mesh',mesh,'figure',[]);
fprintf('PROVISIONAL kayak: mass %.3f kg, draft %.6f m, LWL %.6f m\n', ...
    hs.displacedMass,hs.draft,h.LWL);
disp(h)
try
    [result.resistance,result.details]=delft_resistance(V,h,env,options);
    g=result.details.validity.geometry;
    disp(table(g.names(:),g.values(:),g.lower(:),g.upper(:),g.withinRange(:), ...
        'VariableNames',{'Ratio','Value','Lower','Upper','WithinRange'}))
    result.figure=plot_curve(result,totalMass);
catch exception
    if ~strcmp(exception.identifier,'delft:GeometryRange'), rethrow(exception); end
    result.delftError=exception.message;
    fprintf('Delft rejected this geometry: %s\n',exception.message);
end
end

function fig=plot_curve(r,totalMass)
fig=figure('Name','Provisional kayak resistance','Position',[100 100 1150 550]);
tiledlayout(fig,1,2);
ax=nexttile; moving=r.speed>0;
% Do not draw across the unsupported interval between rest and Fn=0.15.
plot(ax,r.speed(moving),[r.resistance(moving);r.details.Rf(moving);r.details.Rr(moving)].','LineWidth',1.5);
hold(ax,'on'); plot(ax,0,0,'ko','HandleVisibility','off'); yline(ax,0,':','HandleVisibility','off');
xlabel(ax,'Speed through water [m/s]'); ylabel(ax,'Resistance [N]'); grid(ax,'on');
legend(ax,{'Total','Friction','Residuary'},'Location','best');
if r.details.validity.geometryExtrapolation
    titleText='Provisional kayak — Delft geometry extrapolation';
else
    titleText='Provisional kayak — Delft geometry screen passed';
end
sgtitle(fig,titleText);
ax=nexttile; axis(ax,'off'); g=r.details.validity.geometry;
lines={sprintf('Assumed total mass: %.1f kg',totalMass), ...
    'Upright, fixed zero trim; provisional sealed BRI', ...
    'Import conventions unverified; physical accuracy unknown.', ...
    sprintf('Moving speed range: %.3f–%.3f m/s',min(r.speed(moving)),max(r.speed)), ...
    'Zero is isolated; intervening low speeds unsupported.', '', 'Failed geometry checks:'};
for k=find(~g.withinRange)
    lines{end+1}=sprintf('%s: %.6g (limits %.6g–%.6g)',g.names{k},g.values(k),g.lower(k),g.upper(k)); %#ok<AGROW>
end
if g.passed, lines{end+1}='None (screen only, not physical validation).'; end
lines=[lines {'' sprintf('Negative residuary samples: %d',nnz(r.details.validity.negativeResiduary)), ...
    sprintf('Negative total samples: %d',nnz(r.details.validity.negativeTotal))}];
text(ax,0,1,strjoin(lines,newline),'VerticalAlignment','top','Interpreter','none','FontSize',10);
end
