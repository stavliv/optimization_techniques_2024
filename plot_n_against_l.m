function [n_table, lamdas] = plot_n_against_l(functions, a_start, b_start, optim_function)
% plot k against l, fixed e
epsilon = 0.001;
size = 50;

l_end = 100*epsilon;
l_start = 2*epsilon + epsilon/10;
lamdas = linspace(l_start, l_end, size);

n_table = zeros(size, 3);
func_index = 1;
figure;
for func = functions
    % Compute k for all lamdas and store.
    lamdas_index = 1;
    for lamda = lamdas
        if strcmp(optim_function, 'bisection')
            [a_array, b_array, k] = bisection(func, epsilon, lamda, a_start, b_start);
            n = 2 * (k - 1);
        elseif strcmp(optim_function, 'golden_section')
            [a_array, b_array, k] = golden_section(func, lamda, a_start, b_start);
            n = 2 + (k - 1);
        elseif strcmp(optim_function, 'fibonacci_minimize')
            [a_array, b_array, k] = fibonacci_minimize(func, lamda, a_start, b_start);
            n = 2 + (k - 1);
        elseif strcmp(optim_function, 'bisection_derivative')
            [a_array, b_array, k] = bisection_derivative(func, lamda, a_start, b_start);
            n = k - 1;
        else
            error('Invalid optim_function provided, please read the documentation for valid options.');
        end
        n_table(lamdas_index, func_index) = n;
        lamdas_index = lamdas_index + 1;
    end
    % Plot.
    subplot(3, 1, func_index)
    plot(lamdas, n_table(:, func_index), '-r', 'LineWidth', 1.4)
    title(sprintf('$f_%d$', func_index),'Interpreter', 'latex')

    func_index = func_index + 1;
end
xlabel('lamda') 
ylabel('numbers of iterations') 
sgtitle(sprintf('%s', optim_function)) 
