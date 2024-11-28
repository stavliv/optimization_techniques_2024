function [f_min, point_min, points, k] = steepest_descent(f, point_init, epsilon, gamma_mode, gamma_fixed, s_armijo, b_armijo, a_armijo)
% STEEPEST_DESCENT Performs optimization using the steepest descent method.
%
%   This function implements the steepest descent algorithm to find the
%   local minimum of a multivariable function `f`. The method iteratively
%   updates the point in the direction of the negative gradient.
%
%   The step size `gk` can be chosen based on the mode specified in `gamma_mode`.
%   Supported modes:
%     - "fixed": Uses a constant step size `gk`.
%     - "golden_section": Optimizes `gk` using the golden section method.
%     - "armijo": Dynamically selects `gk` using the Armijo rule.
%
% Input Arguments:
%   f - (symbolic function) The multivariable function to minimize.
%   point_init - (1x2 array) Initial point for the optimization.
%   epsilon - (float) Tolerance for the norm of the gradient, the stopping criterion.
%   gamma_mode - (string) Method to determine the step size: "fixed", "golden_section", or "armijo".
%   gamma_fixed - (float) Fixed step size for gradient descent updates (used if `gamma_mode` is "fixed").
%   s_armijo - (float) Initial step size for Armijo rule (required if `gamma_mode` is "armijo").
%   b_armijo - (float) Reduction factor for Armijo rule (required if `gamma_mode` is "armijo").
%   a_armijo - (float) Armijo parameter \( \alpha \) for sufficient decrease (required if `gamma_mode` is "armijo").
%
% Output Arguments:
%   f_min - (float) Minimum value of the function `f` at the final point.
%   point_min - (1x2 array) Coordinates of the final point where `f` is minimized.
%   points - (kx2 matrix) Sequence of points visited during the iterations.
%            Each row represents a point.
%   k - (integer) Number of iterations performed before termination.
%
% Notes:
%   - The function uses symbolic differentiation to compute the gradient.
%   - If `gamma_mode` is "golden_section", the function minimizes `f` along
%     the descent direction at each step to find the best step size `gk` using
%     the golden section method.
%   - If `gamma_mode` is "armijo", it dynamically selects `gk` to ensure sufficient
%     decrease in `f` using the Armijo rule.
%   - The maximum number of iterations is set to 1000.
%
    max_iterations = 1000;

    k = 1;
    point_k = point_init;
    points = [];
    points(1, :) = point_init; % First element assignment in case we don't even get in the loop.
    syms x y
    d = -vpa(subs(jacobian(f), [x, y], point_k));
    d_norm = norm(d);

    % Parameters for the golden section method
    golden_a = 0;  % Start of interval for gk
    golden_b = 2;   % End of interval for gk
    golden_l = epsilon * 2.5;  % Tolerance for golden section search

    while d_norm >= epsilon && k <= max_iterations
        % Compute the gradient and descent direction
        d = -vpa(subs(jacobian(f), [x, y], point_k));
        d_norm = norm(d);

        % Determine the step size based on the gamma mode
        switch gamma_mode
            case "fixed"
                % Use the fixed step size gamma_fixed
                gamma = gamma_fixed;

            case "golden_section"
                % Define the 1D function for line search
                h = @(g) subs(f, [x, y], point_k + g * d);
                gamma = golden_section(h, golden_l, golden_a, golden_b);

            case "armijo"
                % Backtracking line search with Armijo rule
                mk = 0; % Initialize backtracking counter
                gamma = s_armijo * b_armijo^mk; % Initial step size
                while double(subs(f, [x, y], point_k) - subs(f, [x, y], point_k + gamma * d)) < ...
                      double(-a_armijo * gamma * (d' * d))
                    mk = mk + 1; % Increment backtracking step
                    gamma = s_armijo * b_armijo^mk; % Update step size
                end

            otherwise
                error("Unsupported gamma_mode: %s. Supported modes are 'fixed', 'golden_section', and 'armijo'.", gamma_mode);
        end

        % Update the point
        k = k + 1;
        point_k = point_k + gamma * d;
        points(k, :) = point_k;
    end

    % Final results
    point_min = point_k;
    f_min = double(subs(f, [x, y], point_k));
end