% Methods to create Normal Matrix representation for Poission problem
function f = Create_Poisson_problem_A(N)
    A = zeros(N*N, N*N);
    for i = 1: N*N
        row = zeros(1, N*N);
        % each point is system, lets name them vb_2, vb_1, v, vf_1, vf_2
        % vb represents points on the left and vf represents point on the
        % right
        vb_2_index = i - N;
        vb_1_index = i - 1;
        vf_1_index = i + 1;
        vf_2_index = i + N;
    
        if vb_2_index >= 1
           row(vb_2_index) = -1;
        end
        
        skip_prev = mod(i, N) == 1;

        if skip_prev == false && vb_1_index >= 1
            row(vb_1_index) = -1;
        end
    
        row(i) = 4;
    
        skip_next = mod(i, N) == 0;

        if skip_next == false && vf_1_index <= N*N
            row(vf_1_index) = -1;
        end
    
        if vf_2_index <= N*N
            row(vf_2_index) = -1;
        end 
        A(i, :) = row;
    end
    f = A;
end