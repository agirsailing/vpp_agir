function [X,iter,FLAG] = solve_Newton(X0,env,moth,TW, TWS, TWA)
global DATA;

%---------------------------------------------------
% VPP solver unit that solves for VS & TW
%
%
%---------------------------------------------------
step  = ones(5,1)*999;  % Starting values for step vector [VS;TW]
iter  = 0;        % Count iterations
FLAG  = 0;        % Warning flag for maximum iterations
dX =  [1e-4;1e-4;1e-4;1e-4;1e-4]; % Step size for each DOF

X = X0; % The initial guess
% SUGGEST SOLVER LOOP START HERE..
while norm(step)>1e-3

    %-------- KEEP THESE LINES COMMENTED IN EXERCISE 1
    % Count and check number if iterations
    iter = iter+1;          % Count iterations
    if iter>50; FLAG=1; break; end  % Limit for maximum number of iterations
    %-------- KEEP THESE LINES COMMENTED IN EXERCISE 1

    % Calculate 5-DOF residuals vector (5x1)
    R = calc_residuals_Newton(X,env, moth,TW, TWS, TWA);
    
    if any(~isfinite(R)), FLAG = 2; break; end

    fprintf('%2d X=%s  R=%s\n', iter, mat2str(X',4), mat2str(R',3));
    
    % Calculation of derivatives using a finite difference method
    J = zeros(5,5);
    for k = 1:5
        X_pert = X ;
        X_pert(k) = X_pert(k) + dX(k);
        R_pert = calc_residuals_Newton(X_pert, env, moth, TW, TWS, TWA);
        J(:,k) = (R_pert - R) / dX(k);
    end

    if any(~isfinite(J(:))), FLAG = 2; break; end 

    % Newton calculation and limiting
    step =    J \ R;   % Newton step
    max_step = [1.5; 0.05; deg2rad(2); deg2rad(1); deg2rad(1)]; % Limit the step
    step = sign(step) .* min(abs(step),max_step);
    X    = X - step;                        % Next iteration x
    % SUGGEST SOLVER LOOP END HERE
end

% Final check: the residuals themselves must be small, not just the step
if FLAG == 0
    Rf = calc_residuals_Newton(X,env,moth,TW,TWS,TWA);
    if any(abs(Rf(1:4)) > 5) || abs(Rf(5)) > 1e-3
        FLAG = 3;
    end
end