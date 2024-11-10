function [n_table, epsilons] = plot_n_against_e(functions, a_start, b_start, optim_function)


% plot k against e, fixed l
lamda = 0.01;
size = 49;

e_end = 0.0049;
e_start = e_end / size;
epsilons = linspace(e_start, e_end, size);

n_table = zeros(size, length(functions));
func_index = 1;
figure;
for func = functions
    % Compute k for all epsilons and store.
    epsilons_index = 1;
    for epsilon = epsilons
        if optim_function == 'bisection'
            [a_array, b_array, k] = bisection(func, epsilon, lamda, a_start, b_start);
            n = 2 * (k - 1);
        else
            error('Study for epsilons available only for bisection.');
        end
        n_table(epsilons_index, func_index) = n;
        epsilons_index = epsilons_index + 1;
    end
    % Plot.
    subplot(3,1,func_index)
    plot(epsilons, n_table(:, func_index), '-r', 'LineWidth', 1.4)
    title(sprintf('$f_%d$', func_index),'Interpreter', 'latex')

    func_index = func_index + 1;
end
xlabel('\epsilon') 
ylabel('numbers of iterations') 
sgtitle(sprintf('%s', optim_function)) 
