close all; clear all; clc
addpath('utils');
dt = 1;               % Time step
T = 600;

%% Task Initialization

sys.F = 0.95*[1 dt;
              0 1];
sys.G = [1 0;
         0 1];
sys.H = [1 0;
         1 0;
         0 1];
sys.D = 1;
sys.Q = 0.01*[dt^3/3 dt^2/2;
              dt^2/2 dt];
sys.R = [.25    .125   .00125;
         .125   .25    .00125;
         .00125 .00125 .0025];
x0 = [10;0.5];
P0 = [1    1/dt;
      1/dt 2/dt^2];

% Fetch data
[x,y] = getData(x0,T,sys);
n = size(x,1); m = size(y,1); 
mi = 1; p = 3;
% Priors of each sensor
sys1 = sys; sys1.H = sys.H(1,:); sys1.R = sys.R(1,1);
sys2 = sys; sys2.H = sys.H(2,:); sys2.R = sys.R(2,2);
sys3 = sys; sys3.H = sys.H(3,:); sys3.R = sys.R(3,3);

% functions
x_hat_Z = @(Z, n, mu, y) -Z(1:n,n+1:end)*y + mu(1:n) + Z(1:n,n+1:end)*mu(n+1:end);
x_hat_S = @(mu, S, n, y) mu(1:n) + S(1:n, n+1:end) * (S(n+1:end, n+1:end) \ (y - mu(n+1:end)));
P_hat_S = @(S, n) S(1:n, 1:n) - S(1:n, n+1:end) * (S(n+1:end, n+1:end) \ S(n+1:end, 1:n));

%% Sim-No MC

xhat_CI = zeros(n,T);
Phat_CI = zeros(n,n,T);

xhat_MDRO = zeros(n,T);
Phat_MDRO = zeros(n,n,T);

% centralized KF
[xhat_CKF,Phat_CKF] = KF(x0,P0,sys,y);

% CI
x_prev = x0; P_prev = P0;

for i = 1:T
    % sensor 1
    [mu_1, Sigma_1] = predict(x_prev, P_prev, sys1);
    x_1 = x_hat_S(mu_1, Sigma_1, n, y(1,i)); 
    P_1 = P_hat_S(Sigma_1, n);

    % sensor 2
    [mu_2, Sigma_2] = predict(x_prev, P_prev, sys2);
    x_2 = x_hat_S(mu_2, Sigma_2, n, y(2,i)); 
    P_2 = P_hat_S(Sigma_2, n);

    % sensor 3
    [mu_3, Sigma_3] = predict(x_prev, P_prev, sys3);
    x_3 = x_hat_S(mu_3, Sigma_3, n, y(3,i)); 
    P_3 = P_hat_S(Sigma_3, n);

    x_CI = cat(2, x_1, x_2, x_3);
    P_CI = cat(3, P_1, P_2, P_3);
    [xhat_CI(:,i),Phat_CI(:,:,i)] = fusecovint(x_CI,P_CI);
    x_prev = xhat_CI(:,i);
    P_prev = Phat_CI(:,:,i);  
end

%% MDRO
x_prev = x0; P_prev = P0;

method = 'kl_divergence';
parameter = 0.01;

for i = 1:T
    % sensor 1
    [mu_1, Sigma_1] = predict(x_prev, P_prev, sys1);
    % sensor 2
    [mu_2, Sigma_2] = predict(x_prev, P_prev, sys2);
    % sensor 3
    [mu_3, Sigma_3] = predict(x_prev, P_prev, sys3);

    % Concatenate them
    mu = [mu_1, mu_2, mu_3];
    Sigma = cat(3, Sigma_1, Sigma_2, Sigma_3);
    [mu, Sigma] = reconstruct(mu, Sigma, n, mi);

    [S, U, Z, status, objective_value] = MDRO(method, parameter, n, m, mi, p, Sigma, false);
    xhat_MDRO(:,i) = x_hat_Z(Z, n, mu, y(:,i));
    Phat_MDRO(:,:,i) = U;

    x_prev = xhat_MDRO(:,i);
    P_prev = Phat_MDRO(:,:,i);  
end

%%
mse1 = mean(sum((x-xhat_CKF).^2,1));
mse2 = mean(sum((x-xhat_CI).^2,1));
mse3 = mean(sum((x-xhat_MDRO).^2,1));


%% funcs
function [mu_t, Sigma_t] = predict(x_prev, V_prev, sys)
    A_aug = [sys.F; sys.H * sys.F];
    B_aug = [
        sys.G*(sys.Q)*sys.G'                sys.G*(sys.Q)*sys.G'*sys.H'
        sys.H*sys.G*(sys.Q)*sys.G'          sys.H*sys.G*(sys.Q)*sys.G'*sys.H' + sys.R
    ];
    mu_t = A_aug * x_prev;
    Sigma_t = A_aug * V_prev * A_aug' + B_aug;
end