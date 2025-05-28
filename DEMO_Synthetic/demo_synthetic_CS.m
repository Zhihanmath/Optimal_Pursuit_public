% This demo shows how to use OP\CoSaOP\OP-(A)BESS reconstruct the synthetic signal.

close all
clear; clc;

M = 50; % number of measurements
N = 200; % length of the signal
T = 10; % number of nonzero components
test_num = 500; % Independent runs

%% Generate the sparse signal

SNR = 15;         % Signal-to-noise ratio
rng(105,'twister')
q = randperm(N);
beta = zeros(N,1);
beta(q(1:T)) = randn(T,1);

%%

for it = 1:test_num
    fprintf('\n\n Running %d:\n',it);

    % Design Matrix X
    X = rand(M,N);
    X = X./repmat(sqrt((sum(X.^2,1))),[M,1]);
    
    % noiseless signal
    signal = X * beta;
    
    % Generate the noisy observations
    stdnoise = std(signal)*10^(-SNR/20);
    noise = randn(M,1) * stdnoise;
    y = signal + noise;
    
    IND = find(abs(beta)>0);
    
    %% ======================================================================
    %             Algorithm Comparison
    %  ======================================================================
    K = T;
    k_max = 2;
   
    tic
    beta_op(:,it) = op(y, X, K);      
    t_op(it) = toc;

    tic
    beta_op_bess(:,it) = op_bess(y, X, K, k_max);
    t_op_bess(it) = toc;

    tic
    beta_cosaop(:,it) = Cosaop(y, X, K); 
    t_cosaop(it) = toc;
    
    % NMSE
    mse_op(it) = (norm(beta - beta_op(:,it),'fro')/norm(beta,'fro'))^2;
    mse_op_bess(it) = (norm(beta - beta_op_bess(:,it),'fro')/norm(beta,'fro'))^2;
    mse_cosaop(it) = (norm(beta - beta_cosaop(:,it),'fro')/norm(beta,'fro'))^2;
       

    % Ground truth (support of beta)
    support_true = IND;
    non_support_true = find(beta == 0);

    % Support sets for each algorithm
    support_op = find(beta_op(:,it) ~= 0);
    support_op_bess = find(beta_op_bess(:,it) ~= 0);
    support_cosaop = find(beta_cosaop(:,it) ~= 0);

   

    % TPR: Proportion of true positives identified
    tpr_op(it) = length(intersect(support_true, support_op)) / length(support_true);
    tpr_op_bess(it) = length(intersect(support_true, support_op_bess)) / length(support_true);
    tpr_cosaop(it) = length(intersect(support_true, support_cosaop)) / length(support_true);

    % TNR: Proportion of true negatives identified
    tnr_op(it) = length(intersect(non_support_true, setdiff(1:length(beta), support_op))) / length(non_support_true);
    tnr_op_bess(it) = length(intersect(non_support_true, setdiff(1:length(beta), support_op_bess))) / length(non_support_true);
    tnr_cosaop(it) = length(intersect(non_support_true, setdiff(1:length(beta), support_cosaop))) / length(non_support_true);

  
%     % Print
%     fprintf('OP: NMSE: %g,TPR: %g, TNR:%g,time:%g\n',mean(mse_op),mean(tpr_op),mean(tnr_op),mean(t_op));
%     fprintf('OP-(A)BESS: NMSE: %g,TPR: %g, TNR:%g,time:%g\n',mean(mse_op_bess),mean(tpr_op_bess),mean(tnr_op_bess),mean(t_op_bess));
%     fprintf('CoSaOP: NMSE: %g,TPR: %g, TNR:%g,time:%g\n',mean(mse_cosaop),mean(tpr_cosaop),mean(tnr_cosaop),mean(t_cosaop));


end

ind_op=find(tpr_op>=1);
ind_op_bess=find(tpr_op_bess>=1);
ind_cosaop=find(tpr_cosaop>=1);

% Print
fprintf('M=%d, SNR=%d:\n',M,SNR);
fprintf('OP: Successful recovery num:%d, time:%g\n',length(ind_op),mean(t_op));
fprintf('OP-(A)BESS: Successful recovery num:%d,time:%g\n',length(ind_op_bess),mean(t_op_bess));
fprintf('CoSaOP: Successful recovery num:%d,time:%g\n',length(ind_cosaop),mean(t_cosaop));

