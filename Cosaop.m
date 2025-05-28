function x = Cosaop(y, D, K)
    % Cosaop algorithm for sparse signal reconstruction.
    % Parameters:
    %   y - Observation vector (m x 1)
    %   D - Dictionary matrix (m x n)
    %   K - Sparsity level
    % Returns:
    %   x - Sparse coefficient vector (n x 1)
    
    [m, n] = size(D); % Dimensions of the dictionary
    residual = y;     % Initialize residual to the observation vector
    idx_set = [];     % Indices of selected dictionary atoms
    All = 1:n;        % Set of all indices
    x = zeros(n, 1);  % Initialize sparse coefficient vector
    Corr = eye(m);    % Initialize correlation matrix
    Maxiter = 100;    % Maximum number of iterations
    
    for iter = 1:Maxiter
        % 1. Select atoms most correlated with the residual
        projections = zeros(1, n);
        x_old = x;
        for i = 1:n
            projections(i) = (D(:, i)' * residual)^2 / (D(:, i)' * Corr * D(:, i));
        end
        i1 = k_largest_indices(abs(projections), 2 * K);
        union_set = union(idx_set, i1);
        idx_set = union_set;
        res_set = setdiff(All, idx_set);

        % 2. Solve least squares problem with selected atoms
        D_subset = D(:, idx_set);
        cal = inv(D_subset' * D_subset);
        
        % 3. Update sparse coefficients
        coef = zeros(1, length(idx_set));
        for i = 1:length(idx_set)
            [jj_element, sub_matrix, row_vector, col_vector] = process_matrix(cal, i);
            inv_mat = sub_matrix - (col_vector * row_vector) / jj_element;
            idx_set_new = idx_set;
            idx_set_new(i) = [];
            D_subset_new = D(:, idx_set_new);
            x_subset_new = inv_mat * D_subset_new' * y;
            residual_new = y - D_subset_new * x_subset_new;
            coef(i) = norm(residual_new);
        end
        i2 = k_largest_indices(coef, K);
        idx2 = idx_set(i2);
        idx_set = idx2;
        res_set = setdiff(All, idx_set);

        D_subset = D(:, idx_set);
        cal = inv(D_subset' * D_subset);
        x_subset = cal * D_subset' * y;
        Corr = eye(m) - D_subset * cal * D_subset';
        x = zeros(n, 1);
        x(idx_set) = x_subset;

        % 4. Update residual
        residual = y - D_subset * x_subset;

        % 5. Terminate if the residual or x changes are small
        if norm(residual) < 1e-3 || norm(x - x_old) < 1e-3
            %disp(iter);
            break;
        end
    end
end

function indices = k_largest_indices(array, k)
    % Return the indices of the k largest elements in the array
    [~, sorted_indices] = sort(array, 'descend');
    indices = sorted_indices(1:k);
end

function [jj_element, sub_matrix, row_vector, col_vector] = process_matrix(M, j)
    % Process matrix by removing row and column j
    jj_element = M(j, j);
    sub_matrix = M;
    sub_matrix(j, :) = [];
    sub_matrix(:, j) = [];
    row_vector = M(j, :);
    row_vector(j) = [];
    col_vector = M(:, j);
    col_vector(j) = [];
end
