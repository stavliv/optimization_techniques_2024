function [a_array, b_array] = plot_a_b(functions, a_start, b_start, optim_function)

lamdas = [0.0021, 0.01, 0.05, 0.1];

func_index = 1;
for func = functions
    figure;
    lamdas_index = 1;
    for lamda = lamdas
        if strcmp(optim_function, 'bisection')
            [a_array, b_array, k] = bisection(func, 0.001, lamda, a_start, b_start);
        elseif strcmp(optim_function, 'golden_section')
            [a_array, b_array, k] = golden_section(func, lamda, a_start, b_start);
        elseif strcmp(optim_function, 'fibonacci_minimize')
            [a_array, b_array, k] = fibonacci_minimize(func, lamda, a_start, b_start);
        elseif strcmp(optim_function, 'bisection_derivative')
            [a_array, b_array, k] = bisection_derivative(func, lamda, a_start, b_start);
        else
            error('Invalid optim_function provided, please read the documentation for valid options.');
        end
        % Plot.
        subplot(length(lamdas), 1, lamdas_index)
        for i=1:1:k-1
            plot(i, a_array(i), '*r')
            hold on
            plot(i, b_array(i), 'ob')
        end
        xlabel(sprintf('a,b @l=%d', lamda))
        ylabel('[a_k,b_k]')

        lamdas_index = lamdas_index + 1;
    end
    sgtitle(sprintf('%s, %s', optim_function, func)) 

    func_index = func_index + 1;
end
 
