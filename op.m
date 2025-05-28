function x = op(y, D, K)

     % Sparse signal reconstruction using the Orthogonal Matching Pursuit (OP) algorithm.
     %
     % Parameters:
     % y -- Observed signal to be reconstructed, size (m, 1)
     % D -- Dictionary matrix, size (m, n), with columns as dictionary atoms
     % K -- Sparsity level, i.e., the number of non-zero coefficients desired in the reconstructed signal
     %
     % Returns:
     % x -- Sparse coefficient vector, size (n, 1)

    [m, n] = size(D);
    residual = y;          % Initialize residual as the observed signal
    idx_set = [];          % To store the indices of selected dictionary atoms
    x = zeros(n, 1);       % Initialize sparse coefficient vector
    Corr = eye(m);         % Initialize orthogonalization matrix

    for k = 1:K
        % 1. Select the atom that best matches the residual (i.e., the atom with the largest inner product)
        projections = zeros(1, n); 
        for i = 1:n
            numerator = (D(:, i)' * residual)^2;
            denominator = D(:, i)' * Corr * D(:, i);
            projections(i) = numerator / denominator;
        end
        [~, idx] = max(projections); 
        idx_set = [idx_set, idx];    
        
        % 2. Solve the least squares problem over the selected atom subset to update the sparse coefficients
        D_subset = D(:, idx_set);
        cal = inv(D_subset' * D_subset);
        x_subset = cal * (D_subset' * y);
        Corr = eye(m) - D_subset * (cal * D_subset');

        % 3. Update the sparse coefficients
        x(idx_set) = x_subset;

        % 4. Update the residual
        residual = y - D_subset * x_subset;

        % 5. Terminate early if the residual is sufficiently small
        if norm(residual) < 1e-6
            break;
        end
    end
end
