function compare_configurations(nominal, perturbed, names)

params = {'Modulation_efficiency','FWHM','Ratio'};
param_names = ["Modulation efficiency", "\vartheta FWHM", "Nulling ratio"];
nParams = numel(params);
nConfig = numel(names);

for p = 1:nParams
    figure; hold on;

    for i = 1:nConfig

        % Extract the simulation data for the i-th configuration
        simData = perturbed.(params{p})(i, :); 
        
        % Create a boxplot for this configuration; use a small offset for grouping
        boxplot(simData, 'positions', i, 'colors', 'b', 'symbol', 'r+');
        
        % Plot the nominal value as a horizontal line at the configuration's index
        nominalValue = nominal.(params{p})(i);
        plot([i-0.3, i+0.3], [nominalValue, nominalValue], 'r--','LineWidth',1.5);

    end

    set(gca, 'XTick', 1:nConfig, 'XTickLabel', names);
    xlabel('Array Configuration');
    ylabel(param_names{p});
    title(['Distribution of ', params{p}]);
    grid minor;

    if p == 3
        set(gca, "YScale", "log");
    end

    hold off;
    
end



for p = 1:2
    figure; hold on;
    
    violinplot(perturbed.(params{p})');

    for i = 1:nConfig
        % Overlay the nominal value
        nominalValue = nominal.(params{p})(i);
        plot(i, nominalValue, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
    end

    set(gca, 'XTick', categorical(1:nConfig), 'XTickLabel', names);
    xlabel('Array Configuration');

    ylabel(param_names{p});
    title(['Violin plot of ', param_names{p}]);
    grid minor;
    hold off
end


for p = 1:nParams
    figure; hold on;

    for i = 1:nConfig
        simData = perturbed.(params{p})(i, :);  % adjust indexing as needed
        nominalValue = nominal.(params{p})(i);
        % Compute relative differences
        relDiff = (simData - nominalValue) / nominalValue;
        % Scatter plot of simulation index vs. relative difference
        scatter(1:length(simData), relDiff, 'filled');
    end
    xlabel('Simulation Number');
    ylabel(['Relative Difference of ', param_names{p}]);
    title(['Relative Differences for ', param_names{p}]);
    legend(names, 'Location','best');
    grid minor; hold off;
end



for p = 1:nParams
    figure; hold on;

    means = zeros(nConfig, 1);
    stddevs = zeros(nConfig, 1);
    nominalVals = zeros(nConfig, 1);
    for i = 1:nConfig
        simData = perturbed.(params{p})(i, :);
        means(i) = mean(simData);
        stddevs(i) = std(simData);
        nominalVals(i) = nominal.(params{p})(i);
    end
    % Bar plot with error bars
    bar(1:nConfig, means);
    errorbar(1:nConfig, means, stddevs, 'k.', 'LineWidth', 1.5);
    % Overlay nominal values
    plot(1:nConfig, nominalVals, 'ro', 'MarkerSize',8,'LineWidth',2);
    set(gca, 'XTick', 1:nConfig, 'XTickLabel', names);
    xlabel('Array Configuration');
    ylabel(param_names{p});
    title(['Mean and Standard Deviation of ', param_names{p}]);
    legend({'Mean', 'Mean \pm STD','Nominal'}, 'Location','best');

    if p == 3
        set(gca, "YScale", "log")
    end

    grid minor; hold off
end



figure;
for i = 1:nConfig
    for p = 1:nParams
        subplot(nConfig, nParams, (i-1)*nParams + p)
        simData = perturbed.(params{p})(i, :);  % adjust indexing
        plot(simData, 'b.-');
        hold on
        nominalValue = nominal.(params{p})(i);
        yline(nominalValue, 'r--', "LineWidth", 1.5);
        xlabel('Simulation Index');
        ylabel(param_names{p});
        grid minor;
        if i == 1
            title(param_names(p));
        end
        if p == 1
            ylabel([names(i) ' - ' param_names{p}]);
        elseif p == 3
            set(gca, "YScale", "log")
        end
        hold off
    end
end

end