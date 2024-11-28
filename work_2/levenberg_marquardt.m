function [f_min, point_min, points, k] = levenberg_marquardt(f, point_init, epsilon, gamma_mode, gamma_fixed, s_armijo, b_armijo, a_armijo)
% LEVENBERG_MARQUARDT Implements the Levenberg-Marquardt optimization method.
%
%   This function combines gradient descent and Newton's method by modifying
%   the Hessian matrix to improve stability and handle ill-conditioned problems.
%   It minimizes a multivariable function `f`.
%
%   Supported step size modes:
%     - "fixed": Uses a constant step size `gamma`.
%     - "golden_section": Optimizes `gamma` using the golden section method.
%     - "armijo": Dynamically selects `gamma` using the Armijo rule.
%
% Input Arguments:
%   f          - (symbolic function) Function to minimize.
%   point_init - (1x2 array) Initial point for the optimization.
%   epsilon    - (float) Tolerance for the norm of the gradient, the stopping criterion.
%   gamma_mode - (string) Method to determine step size: "fixed", "golden_section", or "armijo".
%   gamma_fixed - (float) Fixed step size (used if `gamma_mode` is "fixed").
%   s_armijo, b_armijo, a_armijo - Parameters for Armijo rule.
%
% Output Arguments:
%   f_min      - (float) Minimum value of the function `f` at the final point.
%   point_min  - (1x2 array) Coordinates of the final point where `f` is minimized.
%   points     - (kx2 matrix) Sequence of points visited during the iterations.
%   k          - (integer) Number of iterations performed before termination.

    max_iterations = 1000;
    k = 1;
    point_k = point_init;
    points = [];
    points(k, :) = point_init; % First element assignment in case we don't even get in the loop.

    syms x y
    gradient_f = jacobian(f, [x, y]); % Symbolic gradient
    hessian_f = hessian(f, [x, y]); % Symbolic Hessian

    g_norm = norm(vpa(subs(gradient_f, [x, y], point_k))); % Initial gradient norm

    lambda_epsilon = 0.1; % Small positive constant added to regularization.

    % Parameters for the golden section method
    golden_a = 0; % Start of interval for gk
    golden_b = 2;  % End of interval for gk
    golden_l = epsilon * 2.5; % Tolerance for the golden section search

    while g_norm >= epsilon && k < max_iterations
        % Compute gradient and Hessian
        grad = vpa(subs(gradient_f, [x, y], point_k));
        hess = vpa(subs(hessian_f, [x, y], point_k));

        % Compute eigenvalues of the Hessian
        eigenvalues = eig(hess);
        max_eig = max(eigenvalues);
        min_eig = min(eigenvalues);

        % Regularize Hessian
        lambda = max_eig + lambda_epsilon; % Start with lambda > max eigenvalue
        updated_hess = hess + lambda * eye(size(hess));

        % Check if positive definite
        if min(eig(updated_hess)) <= 0
            % If not positive definite, use lambda = -min_eig + epsilon, to
            % make for sure positive definite
            lambda = -min_eig + lambda_epsilon;
            updated_hess = hess + lambda * eye(size(hess));
        end

        % Compute descent direction
        dk = -inv(updated_hess) * grad';

        % Determine step size
        switch gamma_mode
            case "fixed"
                gamma = gamma_fixed;

            case "golden_section"
                % Define the 1D function for line search
                h = @(g) subs(f, [x, y], point_k + g * dk');
                gamma = golden_section(h, golden_l, golden_a, golden_b);

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
        points(k, :) = point_k;

        % Update gradient norm
        g_norm = norm(vpa(subs(gradient_f, [x, y], point_k)));
    end

    % Finalize results
    point_min = point_k;
    f_min = double(subs(f, [x, y], point_k));
end