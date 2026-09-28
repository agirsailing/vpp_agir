function [R, details, hull, env] = example_delft_resistance()
%EXAMPLE_DELFT_RESISTANCE Plot documented DSYHS 25 at 10 m waterline length.
% From the repository root: addpath('examples'); example_delft_resistance
[hull, env] = delft_reference_hull();
Fn = linspace(0.15, 0.75, 121);
V = Fn * sqrt(env.g * hull.LWL);
[R, details] = delft_resistance(V, hull, env);
figure('Name', 'Delft 2008: DSYHS 25');
plot(Fn, R, Fn, details.Rf, Fn, details.Rr, 'LineWidth', 1.5);
xlabel('Froude number'); ylabel('Resistance [N]');
legend('Total bare hull', 'Friction (ITTC, LWL)', 'Residuary (Delft 2008)', ...
    'Location', 'northwest');
title('DSYHS 25, LWL = 10 m; reconstructed published geometry');
grid on;
end
