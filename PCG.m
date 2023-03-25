function [x, niters] = PCG(A, b, x0)
% Precondition Matrix A
    L = sparse(ichol(sparse(A), struct('type', 'ict', 'droptol', 1e-3, 'michol', 'off')));
%   M = L * L';
    
    % Initialize x_0 with all zeros
    x = x0;

    r_0 = b - A * x;
    niters = 0;
    r_k = r_0;
    r_k_1 = zeros;
    n_b = norm(b, 2);
    while norm(r_k, 2) > eps(1) * n_b
        %  z = inv(M) * r
        % instead of calculating inverse we'll use traingular solve
        % Mz = r
        % LL' z = r
        % Ly = r
        % L'z = y
        y = L\r_k;
        z_k = L'\y;
        
        if niters == 0
            p_k = r_0;
            p_k_1 = zeros;
        else
            gamma_k = dot(r_k, z_k) / dot(r_k_1, z_k_1);
            p_k = z_k + gamma_k * p_k_1;
        end

        q_k = A * p_k;
        alpha_k = dot(r_k, z_k) / dot(p_k, q_k);
        x = x + alpha_k * p_k;
        r_k_1 = r_k;
        r_k = r_k - alpha_k * q_k;
        z_k_1 = z_k;
        niters = niters + 1;
    end
end