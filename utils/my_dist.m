function [distance] = my_dist(S1, S2, type)
% KL divergence (type 'kl_divergence' or 1) or Wasserstein distance (type 'wasserstein' or 2)
% between N(c,S1) and N(c,S2)

% Map string inputs to corresponding numerical values
if ischar(type) || isstring(type)
    switch lower(type)
        case 'kl_divergence'
            type = 1;
        case 'wasserstein'
            type = 2;
        otherwise
            error('Type undefined.');
    end
end

% Calculate distance based on the type
if type == 1
    distance = trace(S2\S1) - size(S1,1) - log(det(S2\S1));
elseif type == 2
    distance = sqrt(trace(S1+S2-2*sqrtm(sqrtm(S2)*S1*sqrtm(S2))));
else
    error('Type must be either 1 (KL divergence) or 2 (Wasserstein distance).');
end


