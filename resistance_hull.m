function [h, hs] = resistance_hull(mesh, pose, env)
%RESISTANCE_HULL Extract upright static geometry for delft_resistance.
% SI mesh axes: x forward, y starboard, z down. No heel/trim supported.
% Am is immersed section area halfway between waterline endpoints.
% Cp=volume/(LWL*Am); Cm=Am/(BWL*draft). This midship convention is
% explicit; unusual hulls with a different maximum section need review.
% Extraction does not establish applicability to the Delft regression.
if pose.heel ~= 0 || pose.trim ~= 0
    error('hull:ResistancePose','Resistance extraction requires zero heel and trim.');
end
hs=hydrostatics(mesh,pose,env);
if ~strcmp(hs.condition,'partially_submerged') || hs.waterplaneArea<=0
    error('hull:ResistanceWaterplane','A partially submerged hull with a waterplane is required.');
end
v=hs.mesh.vertices; f=hs.mesh.faces;
wp=v(unique(f(hs.mesh.isWaterplane,:)),:);
h=struct('LWL',hs.LWL,'BWL',max(wp(:,2))-min(wp(:,2)), ...
    'draft',hs.draft,'volume',hs.volume,'Swet',hs.wettedArea, ...
    'Awp',hs.waterplaneArea,'xLCB',hs.centreOfBuoyancyWorld(1), ...
    'xLCF',hs.centreOfFlotationWorld(1),'xFPP',hs.waterlineExtent(2));
x=mean(hs.waterlineExtent);
Am=hgeom.immersed_section_area(hs.mesh,x);
if Am<=0, error('hull:ResistanceSection','Midship immersed section has zero area.'); end
h.Cp=h.volume/(h.LWL*Am); h.Cm=Am/(h.BWL*h.draft);
h.midshipArea=Am; h.midshipX=x;
end
