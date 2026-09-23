function [v,f,sections] = section_mesh(sections)
% Shared normalized-perimeter grid, retaining every section vertex.
n=numel(sections); grid=[]; curves=cell(n,1); parameters=cell(n,1);
if n<2 || any(diff([sections.x])<=0)
    error('hull:Sections','At least two strictly increasing stations are required.');
end
for k=1:n
    p=sections(k).yz;
    tol=1e-10*max(max(p)-min(p));
    p=p([true;vecnorm(diff(p),2,2)>tol],:);
    if size(p,1)>1 && norm(p(end,:)-p(1,:))<=tol, p(end,:)=[]; end
    if size(p,1)<3, error('hull:Sections','Degenerate section at station %g.',sections(k).x); end
    % Explicit segment-intersection check; polyshape must not repair input.
    m=size(p,1);
    for i=1:m
        a=p(i,:); b=p(mod(i,m)+1,:);
        for j=i+1:m
            if j==i+1 || (i==1 && j==m), continue; end
            c=p(j,:); d=p(mod(j,m)+1,:);
            if segments_intersect(a,b,c,d,tol)
                error('hull:Sections','Self-intersecting contour at station %g (edges %d and %d).',sections(k).x,i,j);
            end
        end
    end
    area=sum(p(:,1).*p([2:end 1],2)-p([2:end 1],1).*p(:,2))/2;
    if abs(area)<=tol^2, error('hull:Sections','Section has zero area.'); end
    if area<0, p=p([1 end:-1:2],:); end
    % Corresponding contour seam: lowest y, then lowest z.
    [~,order]=sortrows(p,[1 2]); start=order(1);
    p=p([start:end 1:start-1],:);
    sections(k).yz=p;
    p=[p;p(1,:)]; t=[0;cumsum(vecnorm(diff(p),2,2))]; t=t/t(end);
    curves{k}=p; parameters{k}=t; grid=[grid;t]; %#ok<AGROW>
end
% Merge numerically identical fractions without losing meaningful vertices.
grid=uniquetol(grid,1e-12); grid=grid(grid<1-1e-12); m=numel(grid);
v=zeros(n*m,3);
for k=1:n
    p=interp1(parameters{k},curves{k},grid,'linear');
    v((k-1)*m+(1:m),:)=[repmat(sections(k).x,m,1) p];
end
f=zeros(2*m*(n-1),3); row=0;
for k=1:n-1
    for j=1:m
        a=(k-1)*m+j; b=(k-1)*m+mod(j,m)+1; c=a+m; d=b+m;
        row=row+1; f(row,:)=[a b d]; row=row+1; f(row,:)=[a d c];
    end
end
% Constrained triangulation preserves collinear boundary vertices.
for k=[1 n]
    ix=(k-1)*m+(1:m); p=v(ix,2:3);
    dt=delaunayTriangulation(p,[(1:m).' [2:m 1].']);
    if size(dt.Points,1)~=m, error('hull:Sections','Duplicate resampled section points.'); end
    ff=dt.ConnectivityList(isInterior(dt),:);
    for j=1:size(ff,1)
        q=p(ff(j,:),:); a=det([q(2,:)-q(1,:);q(3,:)-q(1,:)]);
        if (k==1 && a>0) || (k==n && a<0), ff(j,[2 3])=ff(j,[3 2]); end
    end
    f=[f;reshape(ix(ff),size(ff))]; %#ok<AGROW>
end
[v,f]=hgeom.validate_mesh(v,f);
end

function yes=segments_intersect(a,b,c,d,tol)
u=b-a; w=d-c;
cross2=@(p,q) p(1)*q(2)-p(2)*q(1);
s=[cross2(u,c-a) cross2(u,d-a)]; t=[cross2(w,a-c) cross2(w,b-c)];
e=tol*max([norm(u) norm(w) tol]);
yes=(s(1)*s(2)<-e^2 && t(1)*t(2)<-e^2);
if yes, return; end
pts=[c;d;a;b]; ends=[a b;a b;c d;c d]; vals=[s t];
for i=1:4
    lo=min(reshape(ends(i,:),2,2).',[],1)-tol;
    hi=max(reshape(ends(i,:),2,2).',[],1)+tol;
    if abs(vals(i))<=e && all(pts(i,:)>=lo) && all(pts(i,:)<=hi), yes=true; return; end
end
end
