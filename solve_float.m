function result = solve_float(hull, loading, options, env)
%SOLVE_FLOAT Solve vertical position or local draft/heel/trim equilibrium.
% options.mode: 'draft' (default) or 'equilibrium'. Angles are degrees.
% loading.mass [kg], loading.cg [x y z] body coordinates [m].
validateattributes(loading.mass,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(loading.cg,{'numeric'},{'size',[1 3],'real','finite'});
if ~isfield(options,'mode'), options.mode='draft'; end
if ~isfield(options,'heel'), options.heel=0; end
if ~isfield(options,'trim'), options.trim=0; end
if ~isfield(options,'angleBounds'), options.angleBounds=[-30 30;-30 30]; end
if ~isfield(options,'maxIterations'), options.maxIterations=30; end
validateattributes(options.heel,{'numeric'},{'scalar','real','finite'});
validateattributes(options.trim,{'numeric'},{'scalar','real','finite'});
if ~isscalar(string(options.mode)), error('hull:Options','Specify one solver mode.'); end
validateattributes(options.angleBounds,{'numeric'},{'size',[2 2],'real','finite'});
validateattributes(options.maxIterations,{'numeric'},{'scalar','integer','nonnegative','finite'});
if any(diff(options.angleBounds,1,2)<=0) || any(abs(options.angleBounds(:))>=89)
    error('hull:Options','Angle bounds must increase and lie strictly inside +/-89 degrees.');
end
if ~ismember(string(options.mode),["draft","equilibrium"])
    error('hull:Options','mode must be draft or equilibrium.');
end
[v,f]=hgeom.validate_mesh(hull.vertices,hull.faces);
hull.vertices=v; hull.faces=f;
scale=max(max(v)-min(v));
q=[options.heel;options.trim];
hgeom.rotation(struct('heel',q(1),'trim',q(2)));
full=hydrostatics(hull,struct('elevation',min(v(:,3))-scale,'heel',0,'trim',0),env);
if loading.mass>=full.displacedMass
    error('hull:Loading','Mass must be below the sealed hull maximum displaced mass (%g kg).',full.displacedMass);
end
if strcmp(options.mode,'equilibrium') && any(q<options.angleBounds(:,1) | q>options.angleBounds(:,2))
    error('hull:Options','Initial angles are outside angleBounds.');
end
step=1e-3; tolerance=1e-8*scale;
[hs,residual,moment]=evaluate(q); iterations=0;
converged=true; message='Vertical force balance converged at prescribed angles.';
if strcmp(options.mode,'equilibrium')
    converged=false; message='Maximum iterations reached.';
    for k=0:options.maxIterations
        iterations=k;
        if norm(residual,inf)<=tolerance
            converged=true; message='Local force and moment equilibrium converged.'; break
        end
        if k==options.maxIterations, break; end
        J=zeros(2);
        for j=1:2
            dq=zeros(2,1); dq(j)=step;
            % One-sided differences near bounds avoid out-of-range evaluations.
            delta=step;
            if q(j)+step>options.angleBounds(j,2), delta=-step; end
            dq(j)=delta;
            if q(j)+delta<options.angleBounds(j,1), delta=(options.angleBounds(j,2)-q(j))/2; dq(j)=delta; end
            [~,r]=evaluate(q+dq); J(:,j)=(r-residual)/delta;
        end
        if rcond(J)<1e-12, message='Singular attitude Jacobian; equilibrium is unresolved.'; break; end
        change=-J\residual; accepted=false;
        for damping=0:20
            candidate=q+change*2^(-damping);
            if any(candidate<options.angleBounds(:,1) | candidate>options.angleBounds(:,2)), continue; end
            [next,r,m]=evaluate(candidate);
            if norm(r)<norm(residual)
                q=candidate; hs=next; residual=r; moment=m; accepted=true; break
            end
        end
        if ~accepted, message='No decreasing step within angle bounds.'; break; end
    end
end
result.pose=hs.pose; result.hydrostatics=hs;
result.converged=converged; result.message=message; result.iterations=iterations;
result.massResidual=hs.displacedMass-loading.mass;
result.forceResidualWorld=[0 0 -env.g*result.massResidual];
result.horizontalResidual=residual.'; result.momentResidualWorld=moment;
result.stability='not_evaluated'; result.restoringStiffness=NaN(2);
result.stiffnessEigenvalues=[NaN;NaN];
if strcmp(options.mode,'equilibrium') && converged
    K=zeros(2);
    for j=1:2
        dq=zeros(2,1); dq(j)=step;
        [~,~,mp]=evaluate(q+dq); [~,~,mm]=evaluate(q-dq);
        tp=generalized(mp,q+dq); tm=generalized(mm,q-dq);
        K(:,j)=-(tp-tm)/(2*deg2rad(step));
    end
    K=(K+K.')/2; values=eig(K);
    result.restoringStiffness=K; result.stiffnessEigenvalues=values;
    threshold=1e-6*loading.mass*env.g*scale;
    if min(values)<-threshold, result.stability='unstable';
    elseif min(values)>threshold, result.stability='stable';
    else, result.stability='near_neutral'; end
end

    function [hs,r,m]=evaluate(angles)
        pose=struct('elevation',0,'heel',angles(1),'trim',angles(2));
        R=hgeom.rotation(pose); w=v*R.';
        bounds=[min(w(:,3))-scale*1e-8 max(w(:,3))+scale*1e-8];
        elevation=fzero(@balance,bounds,optimset('TolX',1e-11*scale,'Display','off'));
        pose.elevation=elevation; hs=hydrostatics(hull,pose,env);
        if abs(hs.displacedMass-loading.mass)>1e-8*loading.mass
            error('hull:Convergence','Draft solve failed its mass-residual tolerance.');
        end
        cg=loading.cg*R.'-[0 0 elevation];
        arm=hs.centreOfBuoyancyWorld-cg; r=arm(1:2).';
        m=cross(arm,[0 0 -hs.buoyancyForce]);
        function value=balance(e)
            pose.elevation=e;
            state=hydrostatics(hull,pose,env);
            value=state.displacedMass-loading.mass;
        end
    end
end

function t=generalized(m,q)
t=[cosd(q(2))*m(1)-sind(q(2))*m(3);m(2)];
end
