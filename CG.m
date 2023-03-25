function [x, niters] = CG(A, b, x0)

    x = x0;
    r = b - A*x;
    p = r;
    niters = 0;
    n_b = norm(b, 2);
    while norm(r, 2) > eps(1) * n_b
        Ap = A*p;
        alpha = (r'*r)/(p'*Ap);
        x = x + alpha*p;
        r_new = r - alpha*Ap;
        beta = (r_new'*r_new)/(r'*r);
        p = r_new + beta*p;
        r = r_new;
        niters = niters + 1;
    end
end