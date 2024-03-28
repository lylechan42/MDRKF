function [S_optimal, U_optimal, Z_optimal, status, objective_value] = DRO(method, parameter, n, m, Sigma, quiet)
    % Initialize parameters
    rho = 0;
    c = 0;
    gamma = 0;
    d = n+m;

    Sigma_half = (sqrtm(Sigma)+sqrtm(Sigma)')/2;

    % Three types of ambiguity sets are allowed
    switch method
        case 'wasserstein'
            rho = parameter;
        case 'kl_divergence'
            c = parameter;
        case 'moment_based'
            gamma = parameter;
        otherwise
            error('Invalid method selected.');
    end
    
    % Set quiet mode for CVX
    if quiet
        cvx_quiet true;
    else
        cvx_quiet false;
    end

    % Common variables and objective
    cvx_begin sdp
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
        subject to
            S >= 0
            S_xx >= 0
            S_yy >= 0
            S - [U, zeros(n,m); zeros(m,n), zeros(m,m)] >= 0 : Z

        % Additional Constraints based on method
        if strcmp(method, 'wasserstein')
            variable V(d, d) symmetric
            trace(S + Sigma - 2*V) <= rho^2
            [Sigma_half*S*Sigma_half, V; V, eye(d)] >= 0
        elseif strcmp(method, 'kl_divergence')
            trace(Sigma\S) - log_det(S) + log_det(Sigma) - d <= c
        elseif strcmp(method, 'moment_based')
            (1-gamma)*Sigma <= S <= (1+gamma)*Sigma
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
