% Sparse Regression
close all
clear; clc;

%% Boston Housing dataset


% Boston
load('data_boston_y.mat')
load('data_boston_X.mat')
y = VarName1;
x = X1;

% % superconductivty :https://archive.ics.uci.edu/dataset/464/superconductivty+data
% load('data_superconductivty_X_y.mat')
% y = train(1:500,end);
% x = train(1:500,1:end-1);

K_big = 5:5:30;


 for it = 1:length(K_big)

    K = K_big(it); % Number of feature
    k_max = 2;
 
    
    tic;
    beta_op(:,it) = op(y, x, K);     
    t_oop(it) = toc;
    
    tic;
    beta_op_bess(:,it) = op_bess(y, x, K, k_max); 
    t_opabess(it) = toc;  
    
    tic;
    beta_cosaop(:,it) = Cosaop(y, x, K); 
    t_cosaop(it) = toc;
    
    
    %% y_hat
      
    y_hat_OP(:,it) = x*beta_op(:,it);
    y_hat_op_bess(:,it) = x*beta_op_bess(:,it);
    y_hat_cosaop(:,it) = x*beta_cosaop(:,it);
    
    %% R^2
    SS_tot = sum((y - mean(y)).^2);  
    
    SS_res_OP(it) = sum((y - y_hat_OP(:,it)).^2);  
    SS_res_op_bess(it) = sum((y - y_hat_op_bess(:,it)).^2); 
    SS_res_cosaop(it) = sum((y - y_hat_cosaop(:,it)).^2);  
    
    R2_OP(it) = 1 - (SS_res_OP(it) / SS_tot);  
    R2_op_bess(it) = 1 - (SS_res_op_bess(it) / SS_tot);  
    R2_cosaop(it) = 1 - (SS_res_cosaop(it)/ SS_tot);  

 end

