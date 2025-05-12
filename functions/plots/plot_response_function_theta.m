function h = plot_response_function_theta(theta_range, RFs, export_setup)
%PLOT_RESPONSE_FUNCTION_THETA Plots the response functions over a range of
%angular coordinates.
%
% INPUTS:
%   theta_range[Np x 1] Vector of angular coordinates [rad]
%   RFs[Ns x Np]        Matrix of response functions, where each row 
%                       corresponds to a different simulation
% OPTIONAL ARGUMENT INPUTS:
%   YScale              Either "linear" or "log" (default: "log").
%   Nominal_curve       Index of the nominal curve (default: 1).
%                      
% OUTPUTS:
%   h[handle]           Handle to the plot
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%   2025-04-10 -------- 1.1
%                     - Added YScale optional argument for plots.
%                     - Result is normalised along the Y axis if desired. 
%                     - Can specify a nominal curve to differentiate style.
%   2025-04-17 -------- 1.2
%                     - Added integrated plot support.
%   2025-05-12 -------- 1.2.1
%                     - Corrected parsing in some case with default
%                       arguments.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    theta_range = default_arguments("theta")
    RFs = []
    export_setup.YScale = "log"
    export_setup.Nominal_curve = 1
    export_setup.Normalize = true
    export_setup.embedded = []
end

if isempty(theta_range)
    theta_range = default_arguments("theta");
end


% Import styling
style_colors;

% Normalisation values
if export_setup.Normalize
    max_val = max(RFs(export_setup.Nominal_curve, :));
else
    max_val = 1;
end

h = figure;
hold on;

% Nominal curve
h1 = plot(rad2mas(theta_range), RFs(export_setup.Nominal_curve, :) / max_val, ...
    "--", "LineWidth", 1.5, "Color", ui_colours(1, :));

% Remove nominal curve from analysis and get colours
RFs(export_setup.Nominal_curve, :) = [];
cols = get_colours(size(RFs, 1), "autumn_night");

for i = 1:size(RFs, 1)
    plot(rad2mas(theta_range), RFs(i, :) / max_val, "LineWidth", 1.5, "Color", cols(i, :));
end

% Create a dummy line for the legend (black, solid)
h_dummy = plot(nan, nan, 'k-', 'LineWidth', 1.5);  

% Leave linear as default
if strcmp(export_setup.YScale, "log")
    set(gca, "YScale", "log");
end

xlabel("Observation angle [mas]")
ylabel("Normalised response function")
legend([h1, h_dummy], {"Nominal", "Perturbed"});
grid minor;

if isstruct(export_setup.embedded)
    export_figures("embedded", export_setup.embedded);
end

end