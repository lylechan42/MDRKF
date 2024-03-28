function isPSD = isPSD(A)
    % Check if matrix A is positive semidefinite

    % Calculate the eigenvalues
    eigenvalues = eig(A);

    % Check if all eigenvalues are greater than or equal to zero
    isPSD = all(eigenvalues >= 0);
end