function [OPD, OPD_max, h] = find_minimum_OPD_single_branch(N, amplitudes, ...
           phases, positions, theta_star, lambda, desired_ratios, export_settings)
%FIND_MINIMUM_OPD_SINGLE_BRANCH Set up an optimisation problem to find, for
%each target nulling depth, the OPD that the second branch should maximally
%display.
%
% INPUTS:
%   N[1]                Number of apertures. [-]
%   amplitudes[Nx1]     Amplitudes associated to each aperture.
%   phases[Nx1]         Vector of phase shifts for each aperture.
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation. [m]
%   desired_ratios[Mx1] Desired nulling ratio(s) at which to compute the
%                       response. If not given, internal values are used.
%   export_setting[struct] Setting for export. 
%
% OUTPUTS:
%   OPD[Mx1]            Vector of OPDs differences for the second branch.
%                       This is the minimum value (see below for 
%                       periodicity).
%   OPD_max[1]          Period.
%   h[1]                Figure handle
%
% NOTES:
%   - OPDs are periodic, therefore each multiple of n * OPD_max, with n
%     integer number will produce the same nulling ratio. 
%   - Perturbations need to be introduced as phase shifts above.
%
% VERSION HISTORY:
%   2025-03-11 -------- 1.0
%   2025-03-13 -------- 2.0
%                     - Function rewritten to use look up tables instead of
%                       numerical solvers. Now works only on single branch.
%   2025-03-27 -------- 2.1
%                     - Adapted to the export settings to save plots.
%   2025-05-12 -------- 2.2
%                     - Correction to the xlabel.
%                     - Added color consistency.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default response: compute the limits below which it the nulling ratio
% cannot go (since it includes the geometric leakage). Assign the interval
% and the number of points to use to generate look-up tables.
[min_ratio, ~] = compute_nulling_ratio(N, amplitudes, phases, ...
                                        positions, theta_star, lambda);
if nargin < 7
    desired_ratios = logspace(log10(min_ratio), 0, 10);
    export_settings = NaN;
    h = NaN;
elseif nargin < 8
    export_settings = NaN;
    h = NaN;
end
OPD = zeros(length(desired_ratios), 1);
points = 1000;

% The phases are bound to 2pi, therefore compute the maximum OPD to have
% the period extension. After OPD_max, everything repeats. Compute the opd
% range of points.
OPD_max = lambda;
opd = logspace(-16, log10(OPD_max), points);

% Generate the map (i.e. compute the nulling ratio for each opd)
for i = 1:points
    opds = [0, opd(i), 0, 0];
    shift_phases = wrapTo2Pi(2*pi*opds/lambda + phases);
    [ratios(i), ~] = compute_nulling_ratio(N, amplitudes, ...
                            shift_phases, positions, theta_star, lambda);
end

% Limit the response to the growing part of OPDs (since in a period, the
% response will grow and decrease). 
[~, i_max] = max(ratios);
opd_range = opd(1:i_max);
ratios_range = ratios(1:i_max);

% Strip unwanted repeating points that block the functioning of interp1:
% extract the indices of unique points, find the counts of each unique
% element, identify non-unique indices and get the indices of non-unique
% elements in the original vector.
[~, ~, idx] = unique(ratios_range, 'stable');
counts = histcounts(idx, 1:max(idx)+1);
non_unique_indices = find(counts > 1);
non_unique_indices_in_ratios = ismember(idx, non_unique_indices);

% Strip the vectors with only non unique elements
ratios_range = ratios_range(~non_unique_indices_in_ratios);
opd_range = opd_range(~non_unique_indices_in_ratios);

% Create the bijective look-up function
opd_fun = @(r) interp1(ratios_range, opd_range, r, "spline", NaN);

% Begin plot if needed
if isstruct(export_settings)
    style_colors;
    h = figure; hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
    xlabel("Nulling ratio")
    ylabel("Minimum OPD [µm]")
    plot(ratios, opd*1e6, "*-", "LineWidth", 1.5, "Color", ui_colours(1, :));
    plot(ratios_range, opd_range*1e6, "-", "LineWidth", 1.5, "Color", ui_colours(9, :));
    grid minor;
end

% For every desired ratio... 
for i = 1:length(desired_ratios)
    
    % Calculate the response...
    desired_ratio = desired_ratios(i);
    req_opd = opd_fun(desired_ratio); % It's periodic! + n * OPD_max;
    OPD(i) = req_opd;
    
    % Verify: compute the value to see the plot result...
    opds = [0, req_opd, 0, 0];
    shift_phases = wrapTo2Pi(2*pi*opds/lambda + phases);
    [rr, ~] = compute_nulling_ratio(N, amplitudes, shift_phases, ...
                                        positions, theta_star, lambda);
    
    % Continue plot
    if isstruct(export_settings)
        xline(desired_ratio, '--', "LineWidth", 1.5, Color=ui_colours(5, :));
    end
end

if isstruct(export_settings)
    export_figures("embedded", export_settings);
end

end
