function [a_array, b_array , k] = golden_section(func, lamda, a_start, b_start)
%GOLDEN_SECTION
    %The interval repeatedly seperated with ratio equal to g=0.618
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
    max_iterations = 1000;
    g = 0.618;
    
    a_array = []; b_array = [];
    a = a_start; b = b_start;
    k = 1;
    a_array(k) = a_start; b_array(k) = b_start;
    
    x1 = a + (1 - g)*(b - a);
    x2 = a + g*(b - a);
    y1 = subs(func, x1);
    y2 = subs(func, x2);
    
    while b - a > lamda
        k = k + 1;
        if y1 > y2
            a = x1;
            x1 = x2;
            y1 = y2; % Avoiding function reevaluation.
            x2 = a + g*(b - a);
            y2 = subs(func, x2);
        else
            b = x2;
            x2 = x1;
            y2 = y1; % Avoiding function reevaluation.
            x1 = a + (1 - g)*(b - a);
            y1 = subs(func, x1);
        end
        a_array(k) = a;
        b_array(k) = b;
    
        if k >= max_iterations
            warning('Maximum number of iterations reached.');
            break;
        end
    end
end