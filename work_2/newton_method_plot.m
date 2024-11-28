clear

syms x y
f(x, y) = x^5 * exp(- x^2 - y^2);

epsilon = 0.001;

% Fixed gamma.
gamma_fixed = 0.4;
% Armijo parameters.
s_armijo = 1;
b_armijo = 0.5;
a_armijo = 0.1;

% Initial points given in the assignment.
% Excluding [-1, 1] since the algorithms get stuck. 
point_inits = [0, 0; 1, -1];
% Gamma modes
gamma_modes = ["fixed", "golden_section", "armijo"];

% Iterate over each initial point
for idx = 1:size(point_inits, 1)
    point_init = point_inits(idx, :);

    % Iterate over each gamma mode
    for gamma_mode = gamma_modes
        % Perform steepest descent with the current gamma mode
        [f_min, point_min, points, k] = ...
            newton_method(f, point_init, epsilon, gamma_mode, gamma_fixed, s_armijo, b_armijo, a_armijo);

        % Plot results for the current run
        plot_results(f, point_init, points, k, gamma_mode);
    end
end
