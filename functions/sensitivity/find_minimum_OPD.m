function delta_OPDs = find_minimum_OPD(N, amplitudes, phases, positions, ...
    theta_star, lambda, autoplot)
%FIND_MINIMUM_OPD Set up an optimisation problem to find, for each target
%nulling depth, the OPDs that each branch should maximally display. 
%
% INPUTS:
%   N[1]                Number of apertures. [-]
%   amplitudes[Nx1]     Amplitudes associated to each aperture.
%   phases[Nx1]         Vector of phase shifts for each aperture.
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation. [m]
%   autoplot[bool]      If given and true, create plots
%
% OUTPUTS:
%   delta_OPDs[mx(N-1)] Vector of OPDs differences for each branch with
%                       respect to the first branch
%
% NOTES:
%   - Perturbations need to be introduced as phase shifts above.
%
% VERSION HISTORY:
%   2025-03-11 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optional inputs
if nargin < 7
    autoplot = false;
end

% Range of desired nulling depths
desired_nulling_depths = logspace(-10, -3, 10); 

% Initialize results
delta_OPDs = zeros(length(desired_nulling_depths), N-1);
mean_delta_OPDs = zeros(length(desired_nulling_depths), 1);

% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'OptimalityTolerance', 1e-20, 'StepTolerance',1e-10);

% Loop through each desired nulling depth
for i = 1:length(desired_nulling_depths)
    
    % Initial guess for relative OPDs (randomized for robustness)
    desired_ratio = desired_nulling_depths(i);
    initial_delta_opds = randn(N-1, 1) * lambda / 10; % Small perturbations

    % Perform optimization
    delta_opds = fmincon(@(delta_opds) objective_function(delta_opds, N, ...
        amplitudes, positions, theta_star, lambda, desired_ratio, phases), ...
        initial_delta_opds, [], [], [], [], -lambda*ones(N-1, 1) * 1e5, lambda*ones(N-1, 1) * 1e5, [], options);
    
    desired_ratio
    compute_nulling_ratio(N, amplitudes, phases, positions, theta_star, lambda)
    compute_nulling_ratio(N, amplitudes, phases + [0, delta_OPDs(1, :)], positions, theta_star, lambda)

    delta_OPDs(i, :) = delta_opds';
    mean_delta_OPDs(i) = mean(delta_opds);
end

% Plot results
if autoplot
    figure;
    for j = 1:N-1
        loglog(desired_nulling_depths, delta_OPDs(:, j), 'DisplayName', ...
            sprintf('\Delta OPD %d', j), 'LineWidth', 1);
        hold on;
    end
    loglog(desired_nulling_depths, mean_delta_OPDs, '--', 'DisplayName', ...
        'Mean \Delta OPD', 'LineWidth', 1.5);
    xlabel('Nulling Depth');
    ylabel('Minimum Relative OPD [m]');
    legend show; grid minor; hold off;
end

end

function error = objective_function(delta_opds, N, amplitudes, positions, ...
                        theta_star, lambda, desired_ratio, ideal_phases)
    
    % Construct full OPDs (first aperture is reference, so OPD_1 = 0)
    opds = [0; delta_opds * 1e-5];

    % Convert OPDs to phases and wrap to 2π
    phases = wrapTo2Pi(2*pi*opds/lambda + ideal_phases(:));

    % Compute nulling ratio
    [ratio, ~] = compute_nulling_ratio(N, amplitudes, phases, ...
                                        positions, theta_star, lambda);

    error = log10(ratio) - log10(desired_ratio);
end
