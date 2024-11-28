clear 
syms x;

f1 = (x-2)^2 + x*log(x+3);
f2 = exp(-2*x) + (x-2)^2;
f3 = exp(x) * (x^3 - 1) + (x - 1)*sin(x);
functions = [f1, f2, f3];

a_start = -1;
b_start = 3;

% Bisection.
% Plot n against e, fixed l.
[n_table, epsilons] = plot_n_against_e(functions, a_start, b_start, 'bisection');
% Plot n against l.
[n_table, lamdas] = plot_n_against_l(functions, a_start, b_start, 'bisection');
% Plot a and b.
[a_array, b_array] = plot_a_b(functions, a_start, b_start, 'bisection');

% Golden section.
% Plot n against l.
[n_table, lamdas] = plot_n_against_l(functions, a_start, b_start, 'golden_section');
% Plot a and b.
[a_array, b_array] = plot_a_b(functions, a_start, b_start, 'golden_section');

% Fibonacci.
% Plot n against l.
[n_table, lamdas] = plot_n_against_l(functions, a_start, b_start, 'fibonacci_minimize');
% Plot a and b.
[a_array, b_array] = plot_a_b(functions, a_start, b_start, 'fibonacci_minimize');

% Bisection derivative.
% Create the derivatives of the original functions
derivative_functions = [diff(f1)];
i = 1;
for func = functions
    derivative_functions(i) = diff(func);
    i = i + 1;
end
% Plot n against l.
[n_table, lamdas] = plot_n_against_l(derivative_functions, a_start, b_start, 'bisection_derivative');
% Plot a and b.
[a_array, b_array] = plot_a_b(derivative_functions, a_start, b_start, 'bisection_derivative');

