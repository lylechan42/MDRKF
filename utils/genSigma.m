function [Sigma, Sigma_star] = genSigma(d)
% genSigma - Generates Sigma and Sigma_star matrices
% Inputs:
%   d - Dimension of the matrices
% Outputs:
%   Sigma - Nominal covariance
%   Sigma_star - Real covariance
%   W(Sigma,Sigma_star) <= d^.5

    % Generate a random matrix A and compute its eigen decomposition
    A = randn(d);
    [R, ~] = eig(A + A');

    % Define lambda and compute Sigma
    lambda = 0.1 + 9.9 * rand(d, 1);
    Sigma = R * diag(lambda) * R';

    Sigma = (Sigma+Sigma')/2; % Ensure Sigma is symmetric

    Sigma_half = R * diag(sqrt(lambda)) * R';

    % Generate another random matrix A_star and compute its eigen decomposition
    A_star = randn(d);
    [R_star, ~] = eig(A_star + A_star');

    % Define lambda_star and compute Delta_star
    lambda_star = rand(d, 1);
    Delta_star = R_star * diag(lambda_star) * R_star';
    Delta_star_half = R_star * diag(sqrt(lambda_star)) * R_star';

    % Compute Sigma_star
    Sigma_star = (Sigma_half + Delta_star_half)^2;
    
    Sigma_star = (Sigma_star+Sigma_star')/2; % Ensure Sigma_star is symmetric

end

