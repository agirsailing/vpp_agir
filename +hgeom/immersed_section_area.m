function Am = immersed_section_area(mesh, x)
%IMMERSED_SECTION_AREA Area of a closed immersed mesh cut at interior x.
% At coincident rings use the aft-side limit; avoid the aft endpoint.
v=mesh.vertices; f=mesh.faces;
integral=0;
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
    q=q-[x mean(v(:,2)) 0];
    integral=integral+q(1,2)*q(2,3)-q(2,2)*q(1,3);
end
Am=abs(integral)/2;
end
