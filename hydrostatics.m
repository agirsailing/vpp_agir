function hs = hydrostatics(hull, pose, env)
%HYDROSTATICS Hydrostatics of a sealed triangular hull at a prescribed pose.
% pose.elevation [m] is upward-positive; heel/trim [deg] use Ry*Rx.
% Body/world axes: x forward, y starboard, z down; world water surface z=0.
% env.water.rho [kg/m^3], env.g [m/s^2]. See HYDROSTATICS.md.
validateattributes(pose.elevation, {'numeric'}, {'scalar','real','finite'});
validateattributes(env.water.rho, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(env.g, {'numeric'}, {'scalar','real','finite','positive'});
[v,f] = hgeom.validate_mesh(hull.vertices,hull.faces);
R = hgeom.rotation(pose);
v = v*R.'-[0 0 pose.elevation];
scale = max(max(v)-min(v)); tol = 1e-10*scale;
v(abs(v(:,3))<=tol,3)=0;
hs.pose=pose; hs.volume=0; hs.displacedMass=0; hs.buoyancyForce=0;
hs.centreOfBuoyancyWorld=[NaN NaN NaN];
hs.centreOfBuoyancyBody=[NaN NaN NaN];
hs.wettedArea=0; hs.waterplaneArea=0;
hs.centreOfFlotationWorld=[NaN NaN 0];
hs.waterplaneMoments=struct('roll',0,'pitch',0,'product',0);
hs.draft=max(0,max(v(:,3))); hs.waterlineExtent=[NaN NaN]; hs.LWL=0;
hs.mesh=struct('vertices',zeros(0,3),'faces',zeros(0,3),'isWaterplane',false(0,1));
if max(v(:,3))<=0
    hs.condition='dry'; return
elseif min(v(:,3))>=0
    hs.condition='fully_submerged';
    sv=v; sf=f; cap=false(size(f,1),1);
else
    hs.condition='partially_submerged';
    % Clip each oriented triangle, keeping z>=0 and its original winding.
    soup=zeros(0,3);
    for k=1:size(f,1)
        p=v(f(k,:),:); q=zeros(0,3);
        if all(p(:,3)==0), continue; end
        for j=1:3
            a=p(j,:); b=p(mod(j,3)+1,:);
            if a(3)>=0, q(end+1,:)=a; end %#ok<AGROW>
            if (a(3)>0 && b(3)<0) || (a(3)<0 && b(3)>0)
                t=a(3)/(a(3)-b(3)); r=a+t*(b-a); r(3)=0;
                q(end+1,:)=r; %#ok<AGROW>
            end
        end
        for j=2:size(q,1)-1
            tri=q([1 j j+1],:);
            if norm(cross(tri(2,:)-tri(1,:),tri(3,:)-tri(1,:)))>tol^2
                soup=[soup;tri]; %#ok<AGROW>
            end
        end
    end
    [sv,~,ix]=uniquetol(soup,tol,'ByRows',true,'DataScale',1);
    sf=reshape(ix,3,[]).'; cap=false(size(sf,1),1);
    % The open edges of the clipped shell trace the waterplane polygons.
    e=[sf(:,[1 2]);sf(:,[2 3]);sf(:,[3 1])];
    [ue,~,g]=unique(sort(e,2),'rows'); counts=accumarray(g,1);
    edges=ue(counts==1,:);
    if any(counts>2) || isempty(edges) || any(abs(sv(unique(edges),3))>tol)
        error('hull:Clip','Could not construct a valid waterplane boundary.');
    end
    if any(accumarray(edges(:),1,[size(sv,1) 1])>2)
        error('hull:Clip','Waterplane touches a non-simple contour; perturb the pose.');
    end
    wp=polyshape();
    while ~isempty(edges)
        loop=edges(1,:); edges(1,:)=[];
        while loop(end)~=loop(1)
            [row,col]=find(edges==loop(end),1);
            if isempty(row), error('hull:Clip','Open waterplane contour.'); end
            loop(end+1)=edges(row,3-col); edges(row,:)=[]; %#ok<AGROW>
        end
        p=sv(loop(1:end-1),1:2);
        wp=xor(wp,polyshape(p(:,1),p(:,2),'Simplify',true,'KeepCollinearPoints',true));
    end
    tr=triangulation(wp); xy=tr.Points; cf=tr.ConnectivityList;
    % Triangles of polyshape triangulation are made outward (towards -z).
    for k=1:size(cf,1)
        p=xy(cf(k,:),:);
        if det([p(2,:)-p(1,:);p(3,:)-p(1,:)])>0
            cf(k,[2 3])=cf(k,[3 2]);
        end
    end
    n=size(sv,1); sv=[sv;xy zeros(size(xy,1),1)];
    sf=[sf;cf+n]; cap=[cap;true(size(cf,1),1)];
    [hs.waterplaneArea,c,I]=waterplane_properties(xy,cf);
    hs.centreOfFlotationWorld=[c 0];
    hs.waterplaneMoments=I;
    hs.waterlineExtent=[min(xy(:,1)) max(xy(:,1))];
    hs.LWL=diff(hs.waterlineExtent);
end
p=sv(sf(:,1),:); q=sv(sf(:,2),:); r=sv(sf(:,3),:);
area=vecnorm(cross(q-p,r-p,2),2,2)/2;
hs.wettedArea=sum(area(~cap));
% Tetrahedra about a nearby origin reduce cancellation after translations.
o=mean(sv,1); a=p-o; b=q-o; c=r-o;
dv=dot(a,cross(b,c,2),2)/6;
hs.volume=sum(dv);
if hs.volume<=0, error('hull:Clip','Submerged mesh has nonpositive volume.'); end
hs.centreOfBuoyancyWorld=o+sum((a+b+c).*dv,1)/(4*hs.volume);
hs.centreOfBuoyancyBody=(hs.centreOfBuoyancyWorld+[0 0 pose.elevation])*R;
hs.displacedMass=env.water.rho*hs.volume;
hs.buoyancyForce=env.g*hs.displacedMass;
hs.mesh=struct('vertices',sv,'faces',sf,'isWaterplane',cap);
end

function [A,c,I]=waterplane_properties(v,f)
a=v(f(:,1),:); b=v(f(:,2),:); d=v(f(:,3),:);
s=abs((b(:,1)-a(:,1)).*(d(:,2)-a(:,2))-(b(:,2)-a(:,2)).*(d(:,1)-a(:,1)))/2;
A=sum(s); c=sum((a+b+d).*s,1)/(3*A);
a=a-c; b=b-c; d=d-c;
xx=sum(s.*(a(:,1).^2+b(:,1).^2+d(:,1).^2+a(:,1).*b(:,1)+a(:,1).*d(:,1)+b(:,1).*d(:,1))/6);
yy=sum(s.*(a(:,2).^2+b(:,2).^2+d(:,2).^2+a(:,2).*b(:,2)+a(:,2).*d(:,2)+b(:,2).*d(:,2))/6);
xy=sum(s.*((a(:,1)+b(:,1)+d(:,1)).*(a(:,2)+b(:,2)+d(:,2))+a(:,1).*a(:,2)+b(:,1).*b(:,2)+d(:,1).*d(:,2))/12);
I=struct('roll',yy,'pitch',xx,'product',xy);
end
