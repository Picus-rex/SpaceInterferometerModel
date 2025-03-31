function perform_statistics(element, analysis, export_setup)
%
% Statistics and Machine Learning Toolbox
arguments
    element (:, :)
    analysis.desired_percentiles = [90, 95, 99]
    export_setup.verbose = true
    export_setup.create_plots = true
    export_setup.embedded = NaN
end

% Np points per simulation, N simulations
[Np, N] = size(element);  

% Preallocate arrays to store metrics
rms_values   = zeros(1, N);
percentiles  = zeros(N, length(analysis.desired_percentiles));      
mu_values    = zeros(1, N);       
sigma_values = zeros(1, N);       

% For every series...
for i = 1:N

    % Extract the data for the current simulation
    data = element(:, i);
    
    % Compute RMS of the OPD for this simulation
    rms_values(i) = sqrt(mean(data.^2));
    
    % Compute the desired percentiles
    percentiles(i, :) = prctile(data, analysis.desired_percentiles);
    
    % Fit a Gaussian (normal distribution) to the data
    % The fitdist function returns a probability distribution object
    pd = fitdist(data, 'Normal');
    mu_values(i)    = pd.mu;
    sigma_values(i) = pd.sigma;
    
    if export_setup.verbose
        fprintf('Simulation %d: RMS = %.4f, 90th = %.4f, 95th = %.4f, 99th = %.4f, mu = %.4f, sigma = %.4f\n', ...
        i, rms_values(i), percentiles(i,1), percentiles(i,2), percentiles(i,3), mu_values(i), sigma_values(i));
    end
end

if export_setup.create_plots

    simIndex = 1:N;
    figure; hold on;
    plot(simIndex, rms_values, 'o-', 'LineWidth', 1.5, 'DisplayName', 'RMS');
    plot(simIndex, mu_values, 's-', 'LineWidth', 1.5, 'DisplayName', 'Mu');
    plot(simIndex, sigma_values, 'd-', 'LineWidth', 1.5, 'DisplayName', 'Sigma');
    xlabel('Simulation Number');
    ylabel('OPD [µm]');
    legend('Location', 'best');
    title('RMS, Mean and Standard Deviation per Simulation');
    grid minor;

    figure; hold on;
    for i = 1:length(analysis.desired_percentiles)
        plot(simIndex, percentiles(:,i), 'o-', 'LineWidth', 1.5, 'DisplayName', string(analysis.desired_percentiles(i)) + 'th Percentile');
    end
    xlabel('Simulation Number');
    ylabel('OPD [µm]');
    legend('Location', 'best');
    title('Percentile Values per Simulation');
    grid minor;

    % Assume RMS is our central measure and we define an "error" as the difference between RMS and the 90th/99th percentiles.
    error_low = rms_values - percentiles(:,1)';
    error_high = percentiles(:,3)' - rms_values;
    
    figure;
    errorbar(simIndex, rms_values, error_low, error_high, 'o-', 'LineWidth', 1.5);
    xlabel('Simulation Number');
    ylabel('OPD [µm]');
    title('RMS with 90th and 99th Percentile Spread');
    grid minor;

    figure;
    scatter(rms_values, sigma_values, 50, 'filled');
    hold on;
    plot([min(rms_values) max(rms_values)], [min(rms_values) max(rms_values)], 'r--', 'LineWidth', 1.5);
    xlabel('RMS [µm]');
    ylabel('Sigma [µm]');
    title('Correlation between RMS and Sigma');
    grid minor;

    figure;
    boxplot(element, 'Labels', arrayfun(@(x) sprintf('Sim %d', x), 1:N, 'UniformOutput', false));
    ylabel('OPD [µm]');
    title('OPD Distribution per Simulation');
    grid minor;

end

end