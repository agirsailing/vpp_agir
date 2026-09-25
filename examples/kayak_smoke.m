% PROVISIONAL geometry smoke check, not verified kayak hydrostatics.
% Assumptions: metres, second header value=x, source z upward, symmetric half hull.
% 'seal' bridges open endpoints, drops interior centreline visits, and orders a
% hull/deck weave as one chain. Original coordinates remain in sourcePoints.
options=struct('lengthScale',1,'xHeaderColumn',2,'axisSigns',[1 1 -1], ...
    'originOffset',[0 0 0],'halfHull',true,'halfHullClosure','seal');
root=fileparts(fileparts(mfilename('fullpath')));
try
    hull=read_bri(fullfile(root,'data','kayak.bri'),options);
catch exception
    if strcmp(exception.identifier,'hull:Sections')
        fprintf('Kayak smoke check rejected unsupported geometry: %s\n',exception.message);
        fprintf('Confirm the export conventions or supply corrected sections before physical use.\n');
        return
    end
    rethrow(exception)
end
env=get_env_params();
z=hull.vertices(:,3);
pose=struct('elevation',(min(z)+max(z))/2,'heel',0,'trim',0);
hs=hydrostatics(hull,pose,env);
fprintf('PROVISIONAL sealed kayak: volume %.6g m^3, wetted area %.6g m^2\n',hs.volume,hs.wettedArea);
fprintf('Sections with closure adjustments: %d / %d\n',sum([hull.sections.closureAdjusted]),numel(hull.sections));
