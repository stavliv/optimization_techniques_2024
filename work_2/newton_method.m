function [f_min, point_min, points, k] = newton_method(f, point_init, epsilon, gamma_mode, gamma_fixed, s_armijo, b_armijo, a_armijo)
% NEWTON_METHOD Performs optimization using the Newton-Raphson method.
%
%   This function implements the Newton-Raphson method to find the
%   local minimum of a multivariable function `f`. It uses the Hessian
%   and gradient of the function to compute the descent direction.
%
%   The step size `gk` can be chosen based on the mode specified in `gamma_mode`.
%   Supported modes:
%     - "fixed": Uses a constant step size `gk`.
%     - "golden_section": Optimizes `gk` using the golden section method.
%     - "armijo": Dynamically selects `gk` using the Armijo rule.
%
% Input Arguments:
%   f          - (symbolic function) The multivariable function to minimize.
%   point_init - (1x2 array) Initial point for the optimization.
%   epsilon    - (float) Tolerance for the norm of the gradient, the stopping criterion.
%   gamma_mode - (string) Method to determine the step size: "fixed", "golden_section", or "armijo".
%   gamma_fixed - (float) Fixed step size for gradient descent updates (used if `gamma_mode` is "fixed").
%   s_armijo, b_armijo, a_armijo - Parameters for Armijo rule.
%
% Output Arguments:
%   f_min      - (float) Minimum value of the function `f` at the final point.
%   point_min  - (1x2 array) Coordinates of the final point where `f` is minimized.
%   points     - (kx2 matrix) Sequence of points visited during the iterations.
%            Each row represents a point.
%   k          - (integer) Number of iterations performed before termination.

    max_iterations = 1000;
    k = 1;
    point_k = point_init;
    points = [];
    points(k, :) = point_init; % First element assignment in case we don't even get in the loop.

    syms x y
    gradient_f = jacobian(f, [x, y]); % Symbolic gradient
    hessian_f = hessian(f, [x, y]); % Symbolic Hessian

    g_norm = norm(vpa(subs(gradient_f, [x, y], point_k)));

    % Parameters for the golden section method
    golden_alpha = 0; % Start of interval for gk
    golden_beta = 2;  % End of interval for gk
    golden_l = epsilon * 2.5; % Tolerance for the golden section search

    while g_norm >= epsilon && k < max_iterations
        % Compute gradient and Hessian
        grad = vpa(subs(gradient_f, [x, y], point_k)); % Gradient
        hess = vpa(subs(hessian_f, [x, y], point_k)); % Hessian matrix

        % Compute the Newton descent direction
        dk = -inv(hess) * grad';

        % Determine step size based on the gamma mode
        switch gamma_mode
            case "fixed"
                gamma = gamma_fixed;

            case "golden_section"
                % Define the 1D function for line search
                h = @(g) subs(f, [x, y], point_k + g * dk');
                gamma = golden_section(h, golden_l, golden_alpha, golden_beta);

            case "armijo"
                % Armijo rule for backtracking line search
                mk = 0; % Backtracking counter
                gamma = s_armijo * b_armijo^mk; % Initial step size
                while double(subs(f, [x, y], point_k) - subs(f, [x, y], point_k + gamma * dk')) < ...
                      double(-a_armijo * gamma * (grad * dk))
                    mk = mk + 1;
                    gamma = s_armijo * b_armijo^mk;
                end

            otherwise
                error("Unsupported gamma_mode: %s. Supported modes are 'fixed', 'golden_section', and 'armijo'.", gamma_mode);
        end

        % Update point
        k = k + 1;
        point_k = point_k + gamma * dk';
        points(k, :) = point_k; % Record the new point

        % Update gradient norm
        g_norm = norm(vpa(subs(gradient_f, [x, y], point_k)));
    end

    % Finalize results
    point_min = point_k;
    f_min = double(subs(f, [x, y], point_k));
end