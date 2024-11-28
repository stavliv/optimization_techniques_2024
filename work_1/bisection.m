function [a_array, b_array, k] = bisection(func, epsilon, lamda, a_start, b_start)
% BISECTION Finds the local minimum of a function within an interval.
    %
    % Description:
    %   This function applies the bisection method to find the local minimum
    %   of the given function `func` within the interval [start_a, start_b].
    %   The method iteratively shrinks the interval based on function values 
    %   at points offset from the midpoint by `epsilon`, until the interval 
    %   length is less than `tolerance`.
    %
    % Input Arguments:
    %   func - Function handle representing the function to minimize.
    %   epsilon - (float) Distance from the midpoint of the interval for evaluation.
    %   lamda - (float) Tolerance for the interval length, the stopping criterion.
    %   a_start - (float) Initial lower bound of the interval.
    %   b_start - (float) Initial upper bound of the interval.
    %
    % Output Arguments:
    %   a_array - Array of lower bounds at each iteration.
    %   b_array - Array of upper bounds at each iteration.
    %   k - (integer) Number of iterations performed. The final interval is
    %   [a_array(k), b_array(k)]
    %
    max_iterations = 1000;
    
    a_array = []; b_array= [];
    a = a_start; b = b_start;
    k = 1;
    a_array(k) = a_start; b_array(k) = b_start;
    
    while b - a > lamda
        mid = (a + b) / 2;
        x1 = mid - epsilon;
        x2 = mid + epsilon;
        k = k + 1;
        if subs(func, x1) < subs(func, x2)
            b = x2;
        else
            a = x1;
        end
        a_array(k) = a;
        b_array(k) = b;
    
        if k >= max_iterations
            warning('Maximum number of iterations reached.');
            break;
        end
    end
end