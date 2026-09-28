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
x=mean(hs.waterlineExtent); integral=0;
% Intersect the closed immersed surface with x=midship. Select the aft
% side at coincident mesh rings so shared edges are counted only once.
for k=1:size(f,1)
    p=v(f(k,:),:);
    if ~(min(p(:,1))<x && max(p(:,1))>=x), continue; end
    q=zeros(0,3);
    for j=1:3
        a=p(j,:); b=p(mod(j,3)+1,:);
        if (a(1)<x && b(1)>=x) || (b(1)<x && a(1)>=x)
            q(end+1,:)=a+(x-a(1))/(b(1)-a(1))*(b-a); %#ok<AGROW>
        end
    end
    if size(q,1)~=2, continue; end
    tangent=cross(cross(p(2,:)-p(1,:),p(3,:)-p(1,:)),[1 0 0]);
    if dot(q(2,:)-q(1,:),tangent)<0, q=q([2 1],:); end
    % Centre coordinates to reduce cancellation for translated meshes.
    q=q-[x mean(wp(:,2)) 0];
    integral=integral+q(1,2)*q(2,3)-q(2,2)*q(1,3);
end
Am=abs(integral)/2;
if Am<=0, error('hull:ResistanceSection','Midship immersed section has zero area.'); end
h.Cp=h.volume/(h.LWL*Am); h.Cm=Am/(h.BWL*h.draft);
h.midshipArea=Am; h.midshipX=x;
end
