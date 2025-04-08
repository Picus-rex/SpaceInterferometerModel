function perform_statistics(element, analysis, export_setup)
%PERFORM_STATISTICS Evaluate the quality of the exported data in terms of
%visualisation of informations on RMS, percentiles and other statistical
%informations on a statistical population.
%
% REQUIRES:
%   Statistics and Machine Learning Toolbox
%
% INPUTS:
%   element[Np x N]     Statistical distribution of Np points over N
%                       series as generated from the plot.
% OPTIONAL ARGUMENTS:
%   desired_percentiles[1xA] Percentiles to evaluate.
%   verbose[bool]       Print results in the command window.
%   create_plots[bool]  Generate plots.
%   embedded[struct]    Global export elements for export_figure() for
%                       every figure (it's a struct array). 
%   label[string]       String of the plotted element.
%   type[string]        Label containing the scale information in format
%                       'prefix_scale', e.g., 'm_1e-6' or 'N_1e9'.
%
% VERSION HISTORY:
%   2025-03-31 -------- 1.0
%   2025-04-01 -------- 1.0.1
%                     - General improvements on the structure of the
%                       function and generalisation of the output.
%   2025-04-03 -------- 1.1
%                     - If type has no prefix (like for the nulling ratio),
%                       then, change a bit the way the units are displayed.
%   2025-04-03 -------- 1.1.1
%                     - Numerically, some elements may be complex,
%                       therefore in plotting only the real part must be
%                       considered.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    element (:, :)
    analysis.desired_percentiles = [90, 95, 99]
    export_setup.verbose = true
    export_setup.create_plots = true
    export_setup.embedded = NaN
    export_setup.label = "OPD"
    export_setup.type = "m_1e-6"
    export_setup.scale = "linear" 
end

% Np points per simulation, N simulations
[Np, N] = size(element);  

% Preallocate arrays to store metrics
rms_values   = zeros(1, N);
percentiles  = zeros(N, length(analysis.desired_percentiles));      
mu_values    = zeros(1, N);       
sigma_values = zeros(1, N);       

% Get visualisation style
[scale, scale_tag] = get_scale_plots(export_setup.type);
if ~strcmp(scale_tag, "")
    elem_label = sprintf("%s [%s]", export_setup.label, scale_tag);
else
    elem_label = sprintf("%s", export_setup.label, scale_tag);
end

% For every series...
for i = 1:N

    % Extract the data for the current simulation
    data = element(:, i);
    
    % Compute RMS of the OPD for this simulation
    rms_values(i) = sqrt(mean(data.^2)) * scale;
    
    % Compute the desired percentiles
    percentiles(i, :) = prctile(data, analysis.desired_percentiles) * scale;
    
    % Fit a Gaussian (normal distribution) to the data
    % The fitdist function returns a probability distribution object
    pd = fitdist(data, 'Normal');
    mu_values(i)    = pd.mu * scale;
    sigma_values(i) = pd.sigma * scale;
    
    if export_setup.verbose
        fprintf('Simulation %d: RMS = %.4f, 90th = %.4f, 95th = %.4f, 99th = %.4f, mu = %.4f, sigma = %.4f\n', ...
        i, rms_values(i), percentiles(i,1), percentiles(i,2), percentiles(i,3), mu_values(i), sigma_values(i));
    end
end

if export_setup.create_plots
    
    % PLOT 1
    simIndex = 1:N;
    
    % Create a grouped bar chart
    figure;
    bar(simIndex, [rms_values', mu_values', sigma_values'], 'grouped');
    
    % Labels & Formatting
    xlabel('Ordered Simulation Index');
    ylabel(elem_label);
    legend({'RMS', 'Fitted mean value', 'Fitted standard deviation'}, 'Location', 'best');
    title('RMS, Mean and Standard Deviation per Simulation');
    grid minor;
    if strcmp(export_setup.scale, "log")
        set(gca, "YScale", "log")
    end

    if isstruct(export_setup.embedded)
        export_figures("embedded", export_setup.embedded(1))
    end
    
    % PLOT 2
    % Aggregate all OPD values across simulations and Compute mean and std
    all_OPD_values = element(:) * scale; 
    global_mu = mean(all_OPD_values);
    global_sigma = std(all_OPD_values);
    
    % Define bins for histogram
    if strcmp(export_setup.scale, "log")
        % Define bin edges in log space
        num_bins = 100;
        edges = log10(logspace(log10(min(all_OPD_values)), log10(max(all_OPD_values)), num_bins));
    else
        % Define bin edges in linear space
        num_bins = 50;
        edges = linspace(min(all_OPD_values), max(all_OPD_values), num_bins);
    end
    
    % Plot histogram of all OPD values
    figure;
    if strcmp(export_setup.scale, "log")
        histogram(real(log10(all_OPD_values)), real(edges), 'Normalization', 'pdf', 'FaceAlpha', 0.6, 'DisplayStyle', 'bar');
    else
        histogram(all_OPD_values, edges, 'Normalization', 'pdf', 'FaceAlpha', 0.6, 'DisplayStyle', 'bar');
    end
    hold on;
    
    % Overlay the global Gaussian distribution
    x_vals = linspace(min(all_OPD_values), max(all_OPD_values), 100);
    y_vals = normpdf(x_vals, global_mu, global_sigma);
    if strcmp(export_setup.scale, "log")
        plot(log10(x_vals), y_vals, 'r-', 'LineWidth', 2, 'DisplayName', 'Global Gaussian Fit');
    else
        plot(x_vals, y_vals, 'r-', 'LineWidth', 2, 'DisplayName', 'Global Gaussian Fit');
    end
    
    % Labels and legend
    xlabel(elem_label);
    ylabel('Probability Density');
    title(sprintf('Aggregated %s Distribution Across Simulations with Gaussian Fit', export_setup.label));
    legend(sprintf('Histogram of %s', export_setup.label), 'Global Gaussian Fit');
    grid on;

    if isstruct(export_setup.embedded)
        export_figures("embedded", export_setup.embedded(2))
    end

    % PLOT 3
    figure; hold on;
    for i = 1:length(analysis.desired_percentiles)
        plot(simIndex, percentiles(:,i), 'o-', 'LineWidth', 1.5, 'DisplayName', string(analysis.desired_percentiles(i)) + 'th Percentile');
    end
    xlabel('Simulation Number');
    ylabel(elem_label);
    legend('Location', 'best');
    title('Percentile Values per Simulation');
    grid minor;
    if strcmp(export_setup.scale, "log")
        set(gca, "YScale", "log")
    end

    if isstruct(export_setup.embedded)
        export_figures("embedded", export_setup.embedded(3))
    end
    
    % PLOT 4
    % Assume RMS is our central measure and we define an "error" as the difference between RMS and the 90th/99th percentiles.
    error_low = rms_values - percentiles(:,1)';
    error_high = percentiles(:,3)' - rms_values;
    
    figure;
    errorbar(simIndex, rms_values, error_low, error_high, 'o-', 'LineWidth', 1.5);
    xlabel('Simulation Number');
    ylabel(elem_label);
    title('RMS with 90th and 99th Percentile Spread');
    grid minor;
    if strcmp(export_setup.scale, "log")
        set(gca, "YScale", "log")
    end

    if isstruct(export_setup.embedded)
        export_figures("embedded", export_setup.embedded(4))
    end
    
    % PLOT 5
    figure;
    scatter(rms_values, sigma_values, 50, 'filled');
    hold on;
    plot([min(rms_values) max(rms_values)], [min(rms_values) max(rms_values)], 'r--', 'LineWidth', 1.5);
    if ~strcmp(scale_tag, "")
        ylabel(sprintf('Standard deviation of fitted distribution [%s]', scale_tag));
        xlabel(sprintf('RMS [%s]', scale_tag));
    else
        ylabel('Standard deviation of fitted distribution');
        xlabel('RMS');
    end
    title('Correlation between RMS and Sigma');
    grid minor;
    if strcmp(export_setup.scale, "log")
        set(gca, "YScale", "log", "XScale", "log")
    end

    if isstruct(export_setup.embedded)
        export_figures("embedded", export_setup.embedded(5))
    end
    
    % PLOT 6
    figure;
    boxplot(element, 'Labels', arrayfun(@(x) sprintf('Sim %d', x), 1:N, 'UniformOutput', false));
    ylabel(elem_label);
    title(sprintf('%s Distribution per Simulation', export_setup.label));
    grid minor;
    if strcmp(export_setup.scale, "log")
        set(gca, "YScale", "log")
    end

    if isstruct(export_setup.embedded)
        export_figures("embedded", export_setup.embedded(6))
    end

end

end