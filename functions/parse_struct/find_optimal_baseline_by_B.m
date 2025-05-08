function B_opt = find_optimal_baseline_by_B(minB, maxB, intensities, ...
    array, apertures_ratio, phase_shifts, theta_star, lambda, autoplot)
%FIND_OPTIMAL_BASELINE_BY_B Find the best optimal baseline to increment the
%nulling ratio.
%
% INPUTS:
%   minB[1]             Minimum baseline to consider. [m]
%   maxB[1]             Maximum baseline to consider. [m]
%   intensities[1xN]    Amplitudes associated to each aperture.
%   array[string]       Type of array (see define_array()).
%   apertures_ratio[1]  Ratio between maximum and minimum size.
%   phase_shifts[1xN]   Phase associated to each aperture [rad]. 
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation (default: 1e-5). [m]
%   autoplot[bool]      If true, visualise the plot (default: true).
%
% OUTPUTS:
%   B_opt[1]            Optimal baseline to increment nulling ratio. [m]
%
% VERSION HISTORY:
%   2025-05-08 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    minB 
    maxB 
    intensities 
    array 
    apertures_ratio 
    phase_shifts 
    theta_star 
    lambda = 1e-5
    autoplot = true
end

% Define and allocate
N = length(phase_shifts);
Bs = linspace(minB, maxB, 1000);
ratios = zeros(length(Bs), 1);

% For every possible baseline...
for i = 1:length(Bs)
    
    % Update the array and compute the nulling ratio
    positions = define_array(array, Bs(i), apertures_ratio);
    ratios(i) = compute_nulling_ratio(N, intensities, phase_shifts, ...
                            positions, theta_star, lambda);
end

if autoplot
    figure;
    semilogy(Bs, ratios, "LineWidth", 1.5)
    xlabel("Maximum baseline [m]")
    ylabel("Nulling ratio")
    grid minor;
end

% Extract the best point.
[~, i_opt] = min(ratios);
B_opt = Bs(i_opt);

end