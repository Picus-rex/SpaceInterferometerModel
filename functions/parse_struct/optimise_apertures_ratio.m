function [optimal_AR, optimal_NR, optimal_IWA] = optimise_apertures_ratio(array, ...
    intensities, phase_shifts, theta_star, lambda, baseline, AR_range, weight)
%OPTIMISE_APERTURES_RATIO Run a solver to extract the best AR in the range
%to minimise both the nulling ratio and the inner working angle.
%
% INPUTS:
%   array[string]       Type of array (see define_array()).
%   intensities[1xN]    Amplitudes associated to each aperture.
%   phase_shifts[1xN]   Phase associated to each aperture [rad]. 
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation. [m]
%   baseline[1]         Maximum baseline of the system. [m]
%   AR_range[2]         Vector with minimum and maximum apertures ratio.
%   weight[1]           Value ∈ [0, 1] to prioritise the nulling ratio over 
%                       the IWA. When it's 0.5, both have the same
%                       importance.
%
% OUTPUTS
%   optimal_AR[1]       Apertures ratio that provides the best metrics in
%                       the given range.
%   optimal_NR[1]       Nulling ratio resulting using the the optimal AR.
%   optimal_IWA[1]      IWA resulting using the optimal AR. [rad]
%
% SEE ALSO:
%   optimise_baseline
%   plot_array_size_sensitivity
%
% VERSION HISTORY:
%   2025-05-15 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimation of values for normalisations
N = length(phase_shifts);
AR_vals = linspace(AR_range(1), AR_range(2), 30);
nulling_vals = zeros(size(AR_vals));
IWA_vals = zeros(size(AR_vals));

for i = 1:length(AR_vals)
    pos = define_array(array, baseline, AR_vals(i));
    nulling_vals(i) = compute_nulling_ratio(N, intensities, phase_shifts, pos, theta_star, lambda);
    [~, baselines] = classify_baselines(intensities, pos, phase_shifts, false);
    IWA_vals(i) = compute_IWA(baselines, lambda);
end

log_nulling = log10(nulling_vals);
log_nulling_norm = @(x) (log10(x) - min(log_nulling)) / (max(log_nulling) - min(log_nulling));
IWA_norm = @(x) (x - min(IWA_vals)) / (max(IWA_vals) - min(IWA_vals));

% Normalised cost function
cost_function = @(AR) objective_function_norm(AR, array, N, intensities, ...
    phase_shifts, theta_star, lambda, baseline, weight, ...
    log_nulling_norm, IWA_norm);
options = optimset('Display', 'off');
[optimal_AR, ~] = fminbnd(cost_function, AR_range(1), AR_range(2), options);

% Extract optimal values
pos = define_array(array, baseline, optimal_AR);
optimal_NR = compute_nulling_ratio(N, intensities, phase_shifts, pos, theta_star, lambda);
[~, baselines] = classify_baselines(intensities, pos, phase_shifts, false);
optimal_IWA = compute_IWA(baselines, lambda);

end





function cost = objective_function_norm(AR, array, N, intensities, phase_shifts, ...
    theta_star, lambda, baseline, weight, log_nulling_norm, IWA_norm)

positions = define_array(array, baseline, AR);
nulling = compute_nulling_ratio(N, intensities, phase_shifts, positions, theta_star, lambda);
[~, baselines] = classify_baselines(intensities, positions, phase_shifts, false);
IWA = compute_IWA(baselines, lambda);

cost = weight * log_nulling_norm(nulling) + (1 - weight) * IWA_norm(IWA);

end
