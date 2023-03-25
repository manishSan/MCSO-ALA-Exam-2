% Method to implement steepest descent algorithm
function [x, niters] = Method_of_Steepest_Descent_ichol(A, b, x0)
    % Precondition Matrix A
    L = sparse(ichol(sparse(A), struct('type', 'ict', 'droptol', 1e-3, ...
        'michol', 'off')));
    % Initialize x_0 with all zeros
    M = L * L';
    cond(M)
    x = x0;
    niters = 0;
    r = b - A * x;
    n_b = norm(b, 2);
    while norm(r, 2) > eps(1) * n_b
       %  p = inv(M) * r
       % instead of calculating inverse we'll use traingular solve
       % Mp = r
       % LL' p = r
       % Lz = r
       % L'p = z
       z = L\r;
       p = L'\z;

       q = A * p;
       alpha = dot(p, r) / dot(p, q);
       x = x + alpha * p;
       r = r - alpha * q;
       niters = niters + 1;
    end
end