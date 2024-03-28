close all; clear all; clc
addpath('utils');
dt = .01;               % Time step
T = 600;
startTime = datetime('now');
% cvx_solver_settings('LOG_DET_IGNORE', 1);

%% Task Initialization

vars.b1 = .5;
vars.b2 = .8;
vars.b3 = .4;
vars.q1 = 1;
vars.q11 = 5;
vars.q12 = 8;
vars.q13 = 10;

[sys,x0,P0] = getSys(dt,3,vars);
coeff = 0.5;
is_TV = true;

% Fetch data
% [x,y] = getRanData(x0,T,sys,coeff,is_TV);
% [x,y] = getData(x0,T,sys);
[x,y] = getCorrData(x0,T,sys,vars);

n = size(x,1); m = size(y,1); 
mi = 1; p = 3;    

% Settings for monte carlo
run_count = 1000;

num_param = 10;
all_c =  3e-1*linspace(0,1,num_param);
all_gamma = 1*linspace(0,1,num_param);


err_CI = zeros(T, run_count);
err_KL_MDRO = zeros(T, run_count, length(all_c));
err_M_MDRO = zeros(T, run_count, length(all_gamma));



% Priors of each sensor
sys1 = sys; sys1.H = sys.H(1,:); sys1.R = sys.R(1,1);
sys2 = sys; sys2.H = sys.H(2,:); sys2.R = sys.R(2,2);
sys3 = sys; sys3.H = sys.H(3,:); sys3.R = sys.R(3,3);

% functions
x_hat_Z = @(Z, n, mu, y) -Z(1:n,n+1:end)*y + mu(1:n) + Z(1:n,n+1:end)*mu(n+1:end);
x_hat_S = @(mu, S, n, y) mu(1:n) + S(1:n, n+1:end) * (S(n+1:end, n+1:end) \ (y - mu(n+1:end)));
P_hat_S = @(S, n) S(1:n, 1:n) - S(1:n, n+1:end) * (S(n+1:end, n+1:end) \ S(n+1:end, 1:n));


%% Sim MC
l1 = length(all_c);
l3 = length(all_gamma);

CI_Time = zeros(run_count,1);
KL_MDRO_Time = zeros(run_count,l1);
M_MDRO_Time = zeros(run_count,l3);
parfor i = 1:run_count
    [temp_CKF,temp_CI] = deal(zeros(T, 1));% for nees of CKF and CI

    [temp_KL1,temp_KL2] = deal(zeros(T, l1));% for mse and nees of KL-MDRO
    [temp_M1,temp_M2] = deal(zeros(T, l3));% for mse and nees of M-MDRO

    temp = zeros(n,T);

    % Get trajectory and observations
    % [x,y] = getRanData(x0,T,sys,coeff,is_TV);
    [x,y] = getCorrData(x0,T,sys,vars);

    % Centralized KF
    [xhat_CKF,Phat_CKF] = KF(x0,P0,sys,y);
    err_CKF(:,i) = sum((x-xhat_CKF).^2,1)';

    % CI
    tic;
    x_prev = x0; P_prev = P0;
    xhat_CI = zeros(n,T);
    Phat_CI = zeros(n,n,T);
    for j = 1:T
        % sensor 1
        [mu_1, Sigma_1] = predict(x_prev, P_prev, sys1);
        x_1 = x_hat_S(mu_1, Sigma_1, n, y(1,j));
        P_1 = P_hat_S(Sigma_1, n);
    
        % sensor 2
        [mu_2, Sigma_2] = predict(x_prev, P_prev, sys2);
        x_2 = x_hat_S(mu_2, Sigma_2, n, y(2,j));
        P_2 = P_hat_S(Sigma_2, n);

        % sensor 3
        [mu_3, Sigma_3] = predict(x_prev, P_prev, sys3);
        x_3 = x_hat_S(mu_3, Sigma_3, n, y(3,j));
        P_3 = P_hat_S(Sigma_3, n);

        x_CI = cat(2, x_1, x_2, x_3);
        P_CI = cat(3, P_1, P_2, P_3);
        [xhat_CI(:,j),Phat_CI(:,:,j)] = fusecovint(x_CI,P_CI);
        x_prev = xhat_CI(:,j);
        P_prev = Phat_CI(:,:,j);
        temp_CI(j) = (x(:,j)-xhat_CI(:,j))'* (Phat_CI(:,:,j) \ (x(:,j)-xhat_CI(:,j)));
    end

    err_CI(:,i) = sum((x-xhat_CI).^2,1)';
    CI_Time(i) = toc;

    % tau-MDRO
    method = 'kl_divergence';
    for a = 1:l1
        tic;
        c = all_c(a);
        x_prev = x0; P_prev = P0;
        xhat = zeros(n,T);
        Phat = zeros(n,n,T);
        for j = 1:T
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

            [S, U, Z, status, objective_value] = MDRO(method, c, n, m, mi, p, Sigma, true);
            xhat(:,j) = x_hat_Z(Z, n, mu, y(:,j));
            Phat(:,:,j) = U;

            x_prev = xhat(:,j);
            P_prev = Phat(:,:,j);
            temp_KL2(j,a) = (x(:,j)-xhat(:,j))'* (Phat(:,:,j) \ (x(:,j)-xhat(:,j)));
        end

        temp_KL1(:,a) = sum((x-xhat).^2,1)';
        KL_MDRO_Time(i,a) = toc;

    end
    err_KL_MDRO(:,i,:) = temp_KL1;


    % M-MDRO
    method = 'moment_based';
    for a = 1:l3
        tic;
        gamma = all_gamma(a);
        x_prev = x0; P_prev = P0;
        xhat = zeros(n,T);
        Phat = zeros(n,n,T);
        for j = 1:T
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

            [S, U, Z, status, objective_value] = MDRO(method, gamma, n, m, mi, p, Sigma, true);
            xhat(:,j) = x_hat_Z(Z, n, mu, y(:,j));
            Phat(:,:,j) = U;

            x_prev = xhat(:,j);
            P_prev = Phat(:,:,j);
            temp_M2(j,a) = (x(:,j)-xhat(:,j))'* (Phat(:,:,j) \ (x(:,j)-xhat(:,j)));
        end

        temp_M1(:,a) = sum((x-xhat).^2,1)';
        M_MDRO_Time(i,a) = toc;
    end
    err_M_MDRO(:,i,:) = temp_M1;

    fprintf('Run i %d is done\n', i);


end

%% Calculate duration
endTime = datetime('now');
executionDuration = endTime - startTime;
disp(['Start time: ', char(startTime)]);
disp(['End time: ', char(endTime)]);
disp(['Execution duration: ', char(executionDuration)]);

%%
% check = squeeze(mean(tmp,3));
CI_Time_avg = mean(CI_Time,1)/T;
KL_MDRO_Time_avg = mean(KL_MDRO_Time,1)/T;

M_MDRO_Time_avg = mean(M_MDRO_Time,1)/T;

mse1 = mean(err_CI,2);
mse3 = mean(err_KL_MDRO,2); [~, minIndex3] = min(mean(mse3,1)); mse_3 = mse3(:,minIndex3);
mse5 = mean(err_M_MDRO,2);  [~, minIndex5] = min(mean(mse5,1)); mse_5 = mse5(:,minIndex5);


%%
% Define a set of colors for the plots
colors = [0, 0.4470, 0.7410;  % Blue
          0.8500, 0.3250, 0.0980;  % Orange
          0.9290, 0.6940, 0.1250;  % Yellow
          0.4940, 0.1840, 0.5560;  % Purple
          0.4660, 0.6740, 0.1880]; % Green

% Plotting
plot(mse1, 'DisplayName', 'Covariance Intersection', 'LineWidth', 2, 'Color', colors(1,:));
hold on;
plot(mse_3, 'DisplayName', 'KL MDRO', 'LineWidth', 2, 'Color', colors(2,:));

plot(mse_5, 'DisplayName', 'Moment-based MDRO', 'LineWidth', 2, 'Color', colors(4,:));
hold off;

% Enhancements
grid on;  % Add grid lines
set(gca, 'FontSize', 12, 'FontName', 'Arial');  % Set font size and type
xlabel('Time/s', 'FontSize', 18);  % X-axis label
ylabel('MSE', 'FontSize', 18);  % Y-axis label
title('Mean Squared Error Comparison', 'FontSize', 16);  % Title
legend('Location', 'best');  % Legend with automatic placement

% Optionally, set axis limits
% Set x-axis ticks and labels for time in seconds
xTicks = 0:100:600;  % Adjust the step and range as needed
xTickLabels = arrayfun(@(x) sprintf('%.1f', x * dt), xTicks, 'UniformOutput', false);
xticks(xTicks);
xticklabels(xTickLabels);

% Optionally, set axis limits
xlim([0, 600]);  % Adjust as needed
ylim([0, max([mse1; mse_3; mse_4; mse_5]) * 1.1]);  % Adjust as needed



%%
plot(all_c,squeeze(mean(mse3,1)));
plot(all_rho,squeeze(mean(mse4,1)));
plot(all_gamma,squeeze(mean(mse5,1)));

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
