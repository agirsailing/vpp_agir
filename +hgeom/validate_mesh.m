function [v,f] = validate_mesh(v,f)
% Require a connected, closed oriented two-manifold; accept global reversal.
validateattributes(v, {'numeric'}, {'2d','ncols',3,'real','finite'});
validateattributes(f, {'numeric'}, {'2d','ncols',3,'integer','positive','<=',size(v,1)});
scale = max(max(v)-min(v));
if size(f,1)<4 || scale<=0
    error('hull:Mesh','Mesh has no enclosed volume.');
end
n = cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
if any(vecnorm(n,2,2)<=eps(scale^2)*32)
    error('hull:Mesh','Mesh contains degenerate triangles.');
end
e = [f(:,[1 2]);f(:,[2 3]);f(:,[3 1])];
[~,~,g] = unique(sort(e,2),'rows');
if any(accumarray(g,1)~=2) || any(accumarray(g,sign(e(:,2)-e(:,1)))~=0)
    error('hull:Mesh','Mesh must be closed with consistent face winding.');
end
G = graph(e(:,1),e(:,2));
if numel(unique(conncomp(G)))~=1 || numel(unique(f(:)))~=size(v,1)
    error('hull:Mesh','Use one connected shell without unused vertices.');
end
w=v-mean(v,1);
vol=sum(dot(w(f(:,1),:),cross(w(f(:,2),:),w(f(:,3),:),2),2))/6;
if abs(vol)<=eps(scale^3)*64
    error('hull:Mesh','Mesh has zero signed volume.');
end
if vol<0, f=f(:,[1 3 2]); end
end
