function [a_array, b_array, k] = bisection_derivative(df, lamda, a_start, b_start)
%BISECTION_DERIVATIVE
    % Input Arguments:
    %   df - Function handle representing the derivative of the function to minimize.
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
    a_array = []; b_array = [];
    a = a_start; b = b_start;
    k = 1;
    a_array(k) = a_start; b_array(k) = b_start;
    
    N = 0; 
    while( N < log2((b - a) / lamda))
        N = N + 1;
    end
    
    for k = 2:(N + 1)
        x1 = (a + b) / 2;
        val = subs(df, x1);
        if val == 0
            return
        elseif val > 0
            b = x1;
        elseif val < 0
            a = x1;
        end
        a_array(k) = a;
        b_array(k) = b;
    end

end