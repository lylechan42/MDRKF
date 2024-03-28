function save_estimations(method, n_sample, alpha, MSE_MDRO_Z, MSE_CI)
    % Save MSE estimations based on the method, sample size, and alpha value

    % Format alpha to match the naming convention (e.g., 0.5 -> 050)
    formatted_alpha = sprintf('%03d', round(alpha * 100));

    % Determine the method prefix for the variable names
    if strcmp(method, 'kl_divergence')
        method_prefix = 'kl';
    elseif strcmp(method, 'wasserstein')
        method_prefix = 'w';
    else
        error('Invalid method specified.');
    end

    % Construct the variable names
    var_name_MDRO = sprintf('%s_%02d_%s_1000_MDRO', method_prefix, n_sample, formatted_alpha);
    var_name_CI = sprintf('%s_%02d_%s_1000_CI', method_prefix, n_sample, formatted_alpha);

    % Create a temporary structure to hold the variables
    tempStruct = struct();

    % Assign the input MSE values to fields in the temporary structure
    tempStruct.(var_name_MDRO) = MSE_MDRO_Z;
    tempStruct.(var_name_CI) = MSE_CI;

    % Check if static_estimation.mat exists, if not, create an empty file
    filename = 'static_estimation.mat';
    if ~exist(filename, 'file')
        save(filename, '-struct', 'tempStruct'); % Creates an empty .mat file with structure
    else
        save(filename, '-struct', 'tempStruct', '-append'); % Append to existing file
    end
end
