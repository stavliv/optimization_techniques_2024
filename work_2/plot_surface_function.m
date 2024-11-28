clear

syms x y
f(x, y) = x^5 * exp(- x^2 - y^2);

figure("Name", sprintf("x^5 * exp(- x^2 - y^2)"));
fsurf(f)
xlabel("x")
ylabel("y")
zlabel("f(x, y)")
title(sprintf("x^5 * exp(- x^2 - y^2)"))

% Find min and max of f and print

% Find the gradient of the function
gradient_f = gradient(f, [x, y]);

% Solve the gradient to find critical points
critical_points = solve(gradient_f == 0, [x, y]);

% Extract critical points as numeric values
critical_x = double(critical_points.x);
critical_y = double(critical_points.y);

% Evaluate the function at the critical points
critical_values = arrayfun(@(cx, cy) double(subs(f, [x, y], [cx, cy])), critical_x, critical_y);

% Find the minimum and maximum values and their locations
[f_min, min_index] = min(critical_values);
[f_max, max_index] = max(critical_values);

min_point = [critical_x(min_index), critical_y(min_index)];
max_point = [critical_x(max_index), critical_y(max_index)];

% Print the minimum and maximum points
fprintf("Minimum point: (x, y) = (%.4f, %.4f), f(x, y) = %.4f\n", min_point(1), min_point(2), f_min);
fprintf("Maximum point: (x, y) = (%.4f, %.4f), f(x, y) = %.4f\n", max_point(1), max_point(2), f_max);
