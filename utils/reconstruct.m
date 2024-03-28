function [mu,Sigma] = reconstruct(mu_XY, Sigma_XY, n, mi)
    % Determine the number of partitions, p, from the size of the third dimension
    [~, ~, p] = size(Sigma_XY);

    % Initialize Sigma with zeros
    Sigma = zeros(n + p * mi, n + p * mi);

    % Initialize mu with zeros
    mu = zeros(n + p * mi, 1);

    % Fill in the Sigma matrix and mu vector
    for i = 1:p
        Y_i_start = n + (i - 1) * mi + 1;
        Y_i_end = n + i * mi;

        % Extract the i-th submatrix and mean vector
        Sigma_XY_i = Sigma_XY(:, :, i);
        mu_XY_i = mu_XY(:, i);

        % Place the submatrix Sigma_XY_i into the correct position in Sigma
        Sigma([1:n, Y_i_start:Y_i_end], [1:n, Y_i_start:Y_i_end]) = Sigma_XY_i;

        % Place the mean vector mu_XY_i into the correct position in mu
        mu([1:n, Y_i_start:Y_i_end]) = mu_XY_i;
    end

    Sigma = (Sigma+Sigma')/2;

end
