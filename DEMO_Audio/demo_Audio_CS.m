close all
clear; clc;
rng(1,'twister') 


% Audio3*: 
info =audioinfo('_0bN5mYLXb0.wav');
[audio,Fs] = audioread('_0bN5mYLXb0.wav');


audiolength = 480;
t = 1:1:audiolength;
audio = audio(:,1);
audio = audio(t);

M = 150; 
N = audiolength;
SNR = 25;  

beta = dct(audio);

IND = find(abs(beta)<=1e-1); 
beta(IND)=zeros(length(IND),1);

K = N-length(IND);
test_num = 100;

for it = 1:test_num
    fprintf('\n\n Running %d:\n',it);
    % Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
    X = randn(M,N);
    X = X./(ones(M,1)*sqrt(sum(X.^2)));

    %% Measurements
    % noiseless signal
    signal = X * beta;


    % Observation noise   
    stdnoise = std(signal)*10^(-SNR/20);
    noise = randn(M,1) * stdnoise;
    y = signal + noise;

    %% ======================================================================
    %             Algorithm Comparison
    %  ======================================================================
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


    fprintf('OP: NMSE: %g, time:%g\n',mean(mse_op),mean(t_op));
    fprintf('OP-(A)BESS: NMSE: %g,time:%g\n',mean(mse_op_bess),mean(t_op_bess));
    fprintf('CoSaOP: NMSE: %g, time:%g\n',mean(mse_cosaop),mean(t_cosaop));

end

std_op = std(mse_op);
std_op_bess = std(mse_op_bess);
std_cosaop = std(mse_cosaop);

mean_op = mean(mse_op);
mean_op_bess = mean(mse_op_bess);
mean_cosaop = mean(mse_cosaop);
