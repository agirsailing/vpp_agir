function result = kayak_resistance(totalMass)
%KAYAK_RESISTANCE Provisional BRI -> loaded hydrostatics -> Delft check.
% totalMass includes kayak, paddler and gear [kg], default 130.
% Import conventions are UNVERIFIED; all 16 sections need sealed closure.
if nargin==0, totalMass=130; end
root=fileparts(fileparts(mfilename('fullpath')));
options=struct('lengthScale',1,'xHeaderColumn',2,'axisSigns',[1 1 -1], ...
    'originOffset',[0 0 0],'halfHull',true,'halfHullClosure','seal');
mesh=read_bri(fullfile(root,'data','kayak.bri'),options);
env=get_env_params();
% CG is a placeholder required by the solver. Fixed-attitude draft depends
% only on mass; this does NOT solve pitch/roll balance or claim a measured CG.
loading=struct('mass',totalMass,'cg',mean(mesh.vertices,1));
floating=solve_float(mesh,loading,struct('mode','draft','heel',0,'trim',0),env);
assert(floating.converged,floating.message);
[h,hs]=resistance_hull(mesh,floating.pose,env);
V=[0 .15 .25 .35 .45 .55 .65 .75]*sqrt(env.g*h.LWL);
result=struct('hull',h,'hydrostatics',hs,'floating',floating,'speed',V, ...
    'resistance',[],'details',[],'delftError','');
fprintf('PROVISIONAL kayak: mass %.3f kg, draft %.6f m, LWL %.6f m\n', ...
    hs.displacedMass,hs.draft,h.LWL);
disp(h)
try
    [result.resistance,result.details]=delft_resistance(V,h,env);
    disp(table(V(:),result.resistance(:),'VariableNames',{'Speed_m_s','Resistance_N'}))
catch exception
    if ~strcmp(exception.identifier,'delft:GeometryRange'), rethrow(exception); end
    result.delftError=exception.message;
    fprintf('Delft rejected this geometry: %s\n',exception.message);
end
end
