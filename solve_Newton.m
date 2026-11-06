function [X,iter,FLAG,FM] = solve_Newton(X0,env,moth,TWS,TWA)
% FLAG: 0 ok, 1 max iterations, 2 NaN/Inf, 4 no step reduces the residuals

iter = 0;
FLAG = 0;
dX   = [1e-3; 1e-3; 1e-4; 1e-4; 1e-4; 1e-3; 1e-3];
tolR = [1; 1; 1; 2; 1e-3; 1e-3; 1e-3];
max_step = [0.3; 0.02; deg2rad(1); deg2rad(0.5); deg2rad(0.5); 0.2; 0.1];

X = X0;
[R,~] = calc_residuals_Newton(X, env, moth, TWS, TWA);
if any(~isfinite(R)), FLAG = 2; end

while FLAG == 0
    if all(abs(R) < tolR), break; end          % converged
    iter = iter + 1;
    if iter > 150, FLAG = 1; break; end

    % Finite-difference Jacobian
    J = zeros(7,7);
    for k = 1:7
        X_pert = X;
        X_pert(k) = X_pert(k) + dX(k);
        [R_pert,~] = calc_residuals_Newton(X_pert, env, moth, TWS, TWA);
        J(:,k) = (R_pert - R) / dX(k);
    end
    if any(~isfinite(J(:))), FLAG = 2; break; end

    % Newton step, scaled as a whole so no component exceeds its limit
    step = J \ R;
    s = min(1, min(max_step ./ max(abs(step), eps)));
    step = s*step;

    % Backtracking: halve the step until the scaled residual norm decreases
    f0 = norm(R ./ tolR);
    t = 1; accepted = false;
    for m = 1:8
        Xn = X - t*step;
        % Ensure state stays in physical limits during line search
        Xn(1) = max(Xn(1), 0.5);
        Xn(2) = min(Xn(2), moth.mv_span + 0.155);
        Xn(3) = min(max(Xn(3), deg2rad(-40)), deg2rad(40));
        Xn(4) = min(max(Xn(4), deg2rad(-15)), deg2rad(15));
        Xn(6) = min(max(Xn(6), 0), 2);
        Xn(7) = min(max(Xn(7), moth.wing_x_pos), moth.wing_x_pos + moth.wing_length);
        [Rn,~] = calc_residuals_Newton(Xn, env, moth, TWS, TWA);
        if all(isfinite(Rn)) && norm(Rn ./ tolR) < f0
            accepted = true; break;
        end
        t = t/2;
    end
    if ~accepted, FLAG = 4; break; end

    X = Xn; R = Rn;
    % fprintf('%3d  f=%.3g  X=%s\n', iter, norm(R./tolR), mat2str(X',3));  % uncomment to debug
end

[~,FM] = calc_residuals_Newton(X, env, moth, TWS, TWA);
%fprintf('%3d f=%.3g X=%s R=%s\n', iter, norm(R./tolR), mat2str(X',3), mat2str(R',3));
end