% Clear existing variables and reset YALMIP
clear;  % Clear workspace
addpath('utils');

% Define the dimensions
mi = 2;     % dim(Y_i)
p = 20;     % # of sensors
m = mi*p;   % dim(Y)
n = 4;      % dim(X)
d = n+m;    % dim(Z=[X;Y])


% M_MDRO setting
rho = .5;

mu = zeros(d, 1); % mean of Z
[Sigma, Sigma_star] = genSigma(d);

Sigma = (Sigma + Sigma') / 2;


% get one z
z = mvnrnd(mu, Sigma_star, 1)';
x_true = z(1:n); % True x
y = z(n+1:end); % One random observation of y

x_hat_Z = @(Z, n, mu, y) -Z(1:n,n+1:end)*y + mu(1:n) + Z(1:n,n+1:end)*mu(n+1:end);
x_hat_S = @(mu, S, n, y) mu(1:n) + S(1:n, n+1:end) * (S(n+1:end, n+1:end) \ (y - mu(n+1:end)));


%% DRO
% [S_optimal, U_optimal, Z_optimal, status, objective_value] = DRO('moment_based',rho, n, m, Sigma, false);


%% fuse with MDRO

[S_optimal, U_optimal, Z_optimal, status, objective_value] = MDRO('moment_based', rho, n, m, mi, p, Sigma, false);

% cvx version has Warning: This linear matrix inequality appears to be unsymmetric. This is
% very likely an error that will produce unexpected results. Please check
% the LMI; and, if necessary, re-enter the model. 

x_hat1 = x_hat_Z(Z_optimal, n, mu, y);

x_hat2 = x_hat_S(mu, S_optimal, n, y);


