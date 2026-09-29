function hull = read_bri(filename, options)
%READ_BRI Read the repository's section-block BRI dialect with explicit mapping.
% Required options: lengthScale, xHeaderColumn (2 or 3), axisSigns (1x3),
% originOffset (1x3, metres, added after scaling/signs), halfHull (logical).
% Header line 2 and per-section scalar trailers are retained uninterpreted.
required={'lengthScale','xHeaderColumn','axisSigns','originOffset','halfHull'};
for k=1:numel(required)
    if ~isfield(options,required{k}), error('hull:Options','Specify options.%s explicitly.',required{k}); end
end
validateattributes(options.lengthScale,{'numeric'},{'scalar','positive','finite','real'});
validateattributes(options.xHeaderColumn,{'numeric'},{'scalar','integer','>=',2,'<=',3});
validateattributes(options.axisSigns,{'numeric'},{'size',[1 3],'real','finite'});
if any(abs(options.axisSigns)~=1), error('hull:Options','axisSigns must contain only +1 or -1.'); end
validateattributes(options.originOffset,{'numeric'},{'size',[1 3],'real','finite'});
validateattributes(options.halfHull,{'logical'},{'scalar'});
if ~isfield(options,'halfHullClosure'), options.halfHullClosure='strict'; end
if ~ismember(string(options.halfHullClosure),["strict","seal"])
    error('hull:Options','halfHullClosure must be strict or seal.');
end
lines=splitlines(string(fileread(filename))); lines=strip(lines);
lines=lines(strlength(lines)>0); pos=1;
if numel(lines)<3, error('hull:BRI','Incomplete BRI header.'); end
hull.name=char(lines(pos)); pos=pos+1;
fileHeader=numeric_line(lines(pos),1,pos); pos=pos+1;
sections=struct('x',{},'yz',{},'sourceHeader',{},'sourceTrailer',{}, ...
    'sourcePoints',{},'closureAdjusted',{}); ended=false;
while pos<=numel(lines)
    head=numeric_line(lines(pos),3,pos); pos=pos+1;
    if all(head==0), ended=true; break; end
    count=head(1);
    if count<3 || count~=fix(count) || pos+count>numel(lines)
        error('hull:BRI','Invalid point count or truncated block at line %d.',pos-1);
    end
    p=zeros(count,2);
    for j=1:count, p(j,:)=numeric_line(lines(pos),2,pos); pos=pos+1; end
    trailer=numeric_line(lines(pos),1,pos); pos=pos+1;
    sourcePoints=p; adjusted=false;
    if options.halfHull
        t=1e-9*max(1,max(abs(p(:))));
        if any(p(:,1)>t) && any(p(:,1)<-t)
            error('hull:BRI','Half-section must stay on one side of the centreline.');
        end
        offCentre=abs(p(1,1))>t || abs(p(end,1))>t;
        if strcmp(options.halfHullClosure,'strict') && offCentre
            error('hull:BRI','Open half-section: confirm geometry or explicitly set halfHullClosure=''seal''.');
        end
        p(abs(p(:,1))<=t,1)=0;
        % Remove zero-area centreline tails only when explicitly sealing.
        if strcmp(options.halfHullClosure,'seal')
            while size(p,1)>2 && all(p(1:2,1)==0)
                adjusted=adjusted || norm(p(1,:)-p(2,:))>t; p(1,:)=[];
            end
            [p,reordered]=separate_hull_deck(p);
            while size(p,1)>2 && all(p(end-1:end,1)==0)
                adjusted=adjusted || norm(p(end,:)-p(end-1,:))>t; p(end,:)=[];
            end
            adjusted=adjusted || offCentre || reordered;
        end
        p=[p;[-p(end:-1:1,1) p(end:-1:1,2)]];
    end
    xyz=[repmat(head(options.xHeaderColumn),size(p,1),1) p];
    xyz=xyz.*(options.lengthScale*options.axisSigns)+options.originOffset;
    sections(end+1)=struct('x',xyz(1,1),'yz',xyz(:,2:3), ...
        'sourceHeader',head,'sourceTrailer',trailer, ...
        'sourcePoints',sourcePoints,'closureAdjusted',adjusted); %#ok<AGROW>
end
if ~ended || pos<=numel(lines), error('hull:BRI','Missing final 0 0 0 marker or unexpected trailing data.'); end
if numel(sections)<2, error('hull:BRI','At least two sections are required.'); end
x=[sections.x];
if all(diff(x)<0), sections=sections(end:-1:1);
elseif any(diff(x)<=0), error('hull:BRI','Stations must be strictly monotonic.'); end
[hull.vertices,hull.faces,hull.sections]=hgeom.section_mesh(sections);
hull.units='m'; hull.axes='x forward, y starboard, z down';
hull.source=struct('filename',char(filename),'header',fileHeader,'options',options);
end

function [p,reordered]=separate_hull_deck(p)
% Order a sealed half-section as hull out to the sheer, then deck back to y=0.
reordered=false;
if size(p,1)<4, return; end
step=vecnorm(diff(p),2,2);
med=median(step); if med==0, med=1e-12; end
thr=5*med; cut=0; ymax=0;
for i=1:numel(step)
    ymax=max(ymax,p(i,1));
    if step(i)>thr && p(i+1,1)<0.25*max(ymax,1e-12), cut=i; break; end
end
if cut==0
    drop=p(:,1)==0; drop([1 end])=false;
    if any(drop), p=p(~drop,:); reordered=true; end
    return
end
hull=p(1:cut,:); tail=p(cut+1:end,:);
deck=tail(tail(:,1)==0,:); rest=tail(tail(:,1)~=0,:);
if isempty(deck), deck=tail(1,:); rest=tail(2:end,:); end
reordered=true;
while ~isempty(rest)
    dh=vecnorm(rest-hull(end,:),2,2); dd=vecnorm(rest-deck(end,:),2,2);
    [gapH,jh]=min(dh); [gapD,jd]=min(dd);
    if norm(hull(end,:)-deck(end,:))<=min(gapH,gapD), break; end
    if gapH<=gapD
        if dd(jh)<gapH
            take=jd; if dh(jd)<dd(jd), take=jh; end
            deck=[deck;rest(take,:)]; rest(take,:)=[]; %#ok<AGROW>
        else
            hull=[hull;rest(jh,:)]; rest(jh,:)=[]; %#ok<AGROW>
        end
    elseif dh(jd)<gapD
        take=jh; if dd(jh)<dh(jh), take=jd; end
        hull=[hull;rest(take,:)]; rest(take,:)=[]; %#ok<AGROW>
    else
        deck=[deck;rest(jd,:)]; rest(jd,:)=[]; %#ok<AGROW>
    end
end
guard=0;
while ~isempty(rest) && guard<size(p,1)+5
    guard=guard+1;
    dh=vecnorm(rest-hull(end,:),2,2); dd=vecnorm(rest-deck(end,:),2,2);
    [gapH,jh]=min(dh); [gapD,jd]=min(dd);
    if min(gapH,gapD)>thr, break; end
    if gapH<=gapD, hull=[hull;rest(jh,:)]; rest(jh,:)=[]; %#ok<AGROW>
    else, deck=[deck;rest(jd,:)]; rest(jd,:)=[]; %#ok<AGROW>
    end
end
if ~isempty(rest)
    error('hull:BRI','Sealed half-section could not be ordered into one contour.');
end
p=[hull;flipud(deck)];
end

function a=numeric_line(s,n,line)
tokens=regexp(char(s),'\s+','split'); a=str2double(tokens);
if numel(a)~=n || any(~isfinite(a))
    error('hull:BRI','Expected %d finite numbers at nonblank line %d.',n,line);
end
end
