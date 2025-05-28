function x = op_bess(y, D, K, k_max)
    % Sparse Signal Reconstruction (ABESS-OP Algorithm)
    % Parameters:
    % y -- Observed signal to be reconstructed, size (m, 1)
    % D -- Dictionary matrix, size (m, n)
    % K -- Sparsity level
    % k_max -- Maximum number of iterations
    
    [m, n] = size(D);
    residual = y;  % Initialize residual as the observed signal
    idx_set = [];  % To store the indices of selected dictionary atoms
    res_set = 1:n; % Remaining atom set
    All = 1:n;
    x = zeros(n, 1); % Initialize sparse coefficient vector

    
    % 1. Initialization
    projections = zeros(1, length(res_set));
    for i = 1:length(res_set)
        projections(i) = (D(:, res_set(i))' * residual) / norm(D(:, res_set(i)));
    end
    idx = k_elements(abs(projections), K);
    idx_set = [idx_set, res_set(idx)];
    res_set(idx) = [];
%      disp('ABESS-OP Initial idx_set:');
%     disp(idx_set);
    
    % 2. Exchange
    while true
        idx_set_new = splice_op(y, D, K, k_max, idx_set, res_set);
%         disp('ABESS-OP New idx_set:');
%         disp(idx_set_new);
        if isequal(idx_set_new, idx_set)
            break;
        else
            idx_set = idx_set_new;
            res_set = setdiff(All, idx_set);
        end
    end

    % 3. Update Coefficients
    D_subset = D(:, idx_set);
    x_subset = (D_subset' * D_subset) \ (D_subset' * y);
    x(idx_set) = x_subset;
end

function idx = k_elements(lst, k)
    [~, sortedIdx] = sort(lst, 'descend');
    idx = sortedIdx(1:k);
end

function idx = get_k_smallest_indices(lst, k)
    [~, sortedIdx] = sort(lst, 'ascend');
    idx = sortedIdx(1:k);
end

function [jj_element, sub_matrix, row_vector, col_vector] = process_matrix(M, j)
    jj_element = M(j, j);
    sub_matrix = M;
    sub_matrix(j, :) = [];
    sub_matrix(:, j) = [];
    row_vector = M(j, :);
    row_vector(j) = [];
    col_vector = M(:, j);
    col_vector(j) = [];
end

function idx_set_new = splice_op(y, D, K, k_max, idx_set, res_set)
    [m, n] = size(D);
    D_subset = D(:, idx_set);
    cal = inv(D_subset' * D_subset);
    x_subset = cal * (D_subset' * y);
    residual = y - D_subset * x_subset;
    Corr = eye(m) - D_subset * cal * D_subset';
    flag = 0;
    
    projections = zeros(1, length(res_set));
    coef = zeros(1, length(idx_set));
    for i = 1:length(res_set)
        projections(i) = (D(:, res_set(i))' * residual)^2 / (D(:, res_set(i))' * Corr * D(:, res_set(i)));
    end
    for i = 1:length(idx_set)
        [jj_element, sub_matrix, row_vector, col_vector] = process_matrix(cal, i);
%         fprintf('sub_matrix %d,%d\n',size(sub_matrix,1),size(sub_matrix,2))
%         fprintf('row_vector %d,%d\n',size(row_vector,1),size(row_vector,2))
%         fprintf('col_vector %d,%d\n',size(col_vector,1),size(col_vector,2))
        inv_mat = sub_matrix - (col_vector * row_vector) / jj_element;
        idx_set_new_tmp = idx_set;
        idx_set_new_tmp(i) = [];
        D_subset_new = D(:, idx_set_new_tmp);
        x_subset_new = inv_mat * (D_subset_new' * y);
        residual_new = y - D_subset_new * x_subset_new;
        coef(i) = norm(residual_new) - norm(residual);
    end
    
    for k = 1:k_max
        i1 = k_elements(abs(projections), k);
        i2 = get_k_smallest_indices(coef, k);
        idx1 = res_set(i1);
        idx2 = idx_set(i2);
        idx_set_new = setdiff(idx_set, idx2);
        idx_set_new = [idx_set_new, idx1];
        
        D_subset_new = D(:, idx_set_new);
        x_subset_new = (D_subset_new' * D_subset_new) \ (D_subset_new' * y);
        residual_new = y - D_subset_new * x_subset_new;

        if norm(residual_new) < norm(residual) - 0.01
            flag = 1;
            break;
        end
    end
    if flag == 0
        idx_set_new = idx_set;
    end
end
