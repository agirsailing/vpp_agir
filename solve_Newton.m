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
    if iter>30; FLAG=1; break; end  % Limit for maximum number of iterations
    %-------- KEEP THESE LINES COMMENTED IN EXERCISE 1

    % Calculate 5-DOF residuals vector (5x1)
    R = calc_residuals_Newton(X,env, moth,TW, TWS, TWA);

    % Calculation of derivatives using a finite difference method
    J = zeros(5,5);
    for k = 1:5
        X_pert = X ;
        X_pert(k) = X_pert(k) + dX(k);
        R_pert = calc_residuals_Newton(X_pert, env, moth, TW, TWS, TWA);
        J(:,k) = (R_pert - R) / dX(k);
    end

    % Newton calculation and limiting
    step =    J \ R;   % Newton step
    max_step = [1.5; 0.05; deg2rad(2); deg2rad(1); deg2rad(1)]; % Limit the step
    step = sign(step) .* min(abs(step),max_step);
    X    = X - step;                        % Next iteration x
    % SUGGEST SOLVER LOOP END HERE
end
