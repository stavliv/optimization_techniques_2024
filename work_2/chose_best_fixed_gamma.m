clear

syms x y
f(x, y) = x^5 * exp(-x^2 - y^2);

epsilon = 0.001;

% Define the range of gamma values to test
gamma_values = 0.1:0.1:1; % Test gamma from 0.1 to 1 in steps of 0.1

% Initial point, we chose [-1, 1]
point_init = [-1, 1];

% Preallocate array for iterations
k_values = zeros(size(gamma_values));

% Test each gamma value
for i = 1:length(gamma_values)
    gamma_fixed = gamma_values(i);

    % Perform steepest descent with the current fixed gamma
    [~, ~, ~, k] = steepest_descent(f, point_init, epsilon, "fixed", gamma_fixed, [], [], []);
    
    % Store the number of iterations
    k_values(i) = k;
end

% Find the gamma that converges the fastest
[min_k, min_k_idx] = min(k_values);
best_gamma = gamma_values(min_k_idx);

% Plot k against gamma
figure("Name", "Iterations vs Fixed Gamma");
plot(gamma_values, k_values, "o-", "LineWidth", 1.5);
xlabel("Fixed Gamma (\gamma)")
ylabel("Iterations (k)")
title("Convergence Speed vs Fixed Gamma")
grid on

% Highlight the best gamma on the plot
hold on
plot(best_gamma, min_k, "ro", "MarkerFaceColor", "red", ...
     "DisplayName", sprintf("Best Gamma: %.2f (k = %d)", best_gamma, min_k));
legend show
hold off

% Display the best gamma in the console
fprintf("Best Gamma: %.2f (Converged in %d iterations)\n", best_gamma, min_k);