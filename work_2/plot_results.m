function plot_results(f, point_init, points, k, gamma_mode)
    % Calculates function values and generates plots for minimization results.
    %
    % Inputs:
    %   f          - symbolic function being minimized
    %   point_init - initial point of the descent
    %   points     - trajectory points from steepest descent
    %   k          - number of iterations
    %   gamma_mode - mode for selecting step size ("fixed", "golden_section", "armijo")
    %
    % Outputs:
    %   Generates two plots:
    %       1. 3D trajectory on the function surface
    %       2. Function values vs iterations plot
    
    % Calculate function values along the trajectory
    f_values = arrayfun(@(i) double(subs(f, [sym('x'), sym('y')], points(i, :))), 1:k);

    % 3D Trajectory Plot

    figure("Name", sprintf("Trajectory - Init (%.1f, %.1f) - Mode: %s", ...
                           point_init(1), point_init(2), gamma_mode));
    fsurf(f, [-2, 2, -2, 2]) % Adjust the range as necessary for better visibility
    xlabel("x")
    ylabel("y")
    zlabel("f(x, y)")
    title(sprintf("Trajectory for Init Point (%.1f, %.1f) - Mode: %s", ...
                  point_init(1), point_init(2), gamma_mode))
    hold on

    % Plot trajectory on the 3D function surface
    plot3(points(:, 1), points(:, 2), f_values, "o-", ...
          "DisplayName", "Path");

    % Mark the initial point
    plot3(point_init(1), point_init(2), double(subs(f, [sym('x'), sym('y')], point_init)), "ro", ...
          "MarkerFaceColor", "red", "DisplayName", "Initial Point");

    legend show
    hold off

    % Function Values vs. Iterations Plot

    figure("Name", sprintf("Function Values vs Iterations - Init (%.1f, %.1f) - Mode: %s", ...
                           point_init(1), point_init(2), gamma_mode));
    plot(1:k, f_values, "o-", "DisplayName", "f(x, y) per Iteration");
    xlabel("Iteration (k)")
    ylabel("f(x, y)")
    title(sprintf("Convergence for Init Point (%.1f, %.1f) - Mode: %s", ...
                  point_init(1), point_init(2), gamma_mode))
    % Add legend with iteration count and final minimum value
    legend(sprintf("Iterations (k): %d\nFinal Minimum Value: %.4f", k, f_values(end)));
    
    legend show
    hold off
end