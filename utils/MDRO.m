function [S_optimal, U_optimal, Z_optimal, status, objective_value] = MDRO(method, parameter, n, m, mi, p, Sigma, quiet)
% Marginal distributionally robust optmization with Wasserstein distance 
% Inputs:
%   method - Wasserstein or KL divergence or Moment-based set
%   parameter - radius of ambiguity set
%   n - dim x
%   m - dim y
%   mi - dim yi
%   p - # of yi(sensors)
%   Sigma - nominal covariance
%   quiet - cvx_quiet
% Outputs:
%   S_optimal - covaraince of P_{x,y} after MDRO
%   U_optimal - S_xx - S_xy S_yy^-1 S_yx, i.e. covariance of estimator
    
    if isscalar(parameter)
        param = repmat(parameter, p, 1); % Replicate the scalar for each sensor
    elseif isvector(parameter) && length(parameter) == p
        param = parameter; % Use the provided list
    else
        error('Parameter must be a scalar or a list of length equal to p.');
    end

    
    % Set quiet mode for CVX
    if quiet
        cvx_quiet true;
    else
        cvx_quiet false;
    end

    % CVX optimization problem
    cvx_begin sdp
        % Define decision variables
        variable U(n, n) symmetric
        variable S_xx(n, n) symmetric
        variable S_yy(m, m) symmetric
        variable S_xy(n, m)
        
        dual variable Z

        % Define the block matrix S
        S = [S_xx, S_xy; S_xy', S_yy];

        % Objective function
        minimize(-trace(U))

        % Common Constraints
        S >= 0
        S_xx >= 0
        S_yy >= 0
        % Schur complement constraint
        S - [U, zeros(n,m); zeros(m,n), zeros(m,m)] >= 0 : Z

        % Additional Constraints based on method

        if strcmp(method, 'kl_divergence')

            % Loop over each sensor for Wasserstein distance constraints
            for i = 1:p
                Y_i_start = n + (i-1)*mi + 1;
                Y_i_end = n + i*mi;
                Sigma_XY_i = Sigma([1:n, Y_i_start:Y_i_end], [1:n, Y_i_start:Y_i_end]);

                % Sensor-specific constraints
                trace(Sigma_XY_i\S([1:n, Y_i_start:Y_i_end], [1:n, Y_i_start:Y_i_end])) - log_det(S([1:n, Y_i_start:Y_i_end], [1:n, Y_i_start:Y_i_end])) + log_det(Sigma_XY_i) - (n+mi) <= param(i)
            end

        elseif strcmp(method, 'moment_based')

            % Loop over each sensor for Wasserstein distance constraints
            for i = 1:p
                Y_i_start = n + (i-1)*mi + 1;
                Y_i_end = n + i*mi;
                Sigma_XY_i = Sigma([1:n, Y_i_start:Y_i_end], [1:n, Y_i_start:Y_i_end]);

                % Sensor-specific constraints
                (1-param(i))*Sigma_XY_i <= S([1:n, Y_i_start:Y_i_end], [1:n, Y_i_start:Y_i_end]) <= (1+param(i))*Sigma_XY_i
            end

        end
        
        

    cvx_end 

    % After solving the optimization problem
    status = cvx_status;
    objective_value = cvx_optval;

    % Extract optimal values
    S_optimal = S;
    U_optimal = U;
    Z_optimal = Z;
end
