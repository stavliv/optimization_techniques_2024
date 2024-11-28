function [a_array, b_array, k] = fibonacci_minimize(func, lamda, a_start, b_start)
%FIBONACCI_MINIMIZE
    %
    % Input Arguments:
    %   func - Function handle representing the function to minimize.
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
    lamda_intervals = (b_start - a_start) / lamda; 
    N = 0;
    while lamda_intervals > fibonacci(N+1)   
        N = N + 1;
    end
    
    a_array = []; b_array = [];
    a = a_start; b = b_start;
    k = 1;
    a_array(k) = a_start; b_array(k) = b_start;
    
    c = a + (1 - fibonacci(N-1-k)/fibonacci(N-k)) * (b - a);
    d = a + (fibonacci(N-1-k)/fibonacci(N-k)) * (b - a);
    yc = subs(func, c);
    yd = subs(func, d);
    
    for k = 2:N-1
        if yc <= yd
            b = d;  
            d = c;
            yd = yc; % Avoiding function reevaluation.
            c =  a + (1 - fibonacci(N-1-k)/fibonacci(N-k) )*(b - a);
            yc = subs(func, c);
        else
            a = c; 
            c = d; 
            yc = yd; % Avoiding function reevaluation.
            d =  a + (fibonacci(N-1-k)/fibonacci(N-k))*(b - a);
            yd = subs(func, d);
        end
        a_array(k) = a;
        b_array(k) = b;
    end
end