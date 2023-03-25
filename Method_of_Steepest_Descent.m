% Method to implement steepest descent algorithm
function [x, niters] = Method_of_Steepest_Descent(A, b, x0)
    % Initialize x_0 with all zeros

    cond(A)

    x = x0;
    niters = 0;
    r = b - A * x;
    n_b = norm(b, 2);
    while norm(r, 2) > eps(1) * n_b
       p = r;
       q = A * p;
       alpha = dot(p, r) / dot(p, q);
       x = x + alpha * p;
       r = r - alpha * q;
       niters = niters + 1;
    end
end
