% Retrieve the covariance matrix for prespecified position
function[sigma] = GetSigmaForRegime(input_sigma, T, position)

sigma = input_sigma (1 + T(1,2)*(position-1) : T(1,2) * position, :) ;

