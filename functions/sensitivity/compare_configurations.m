function compare_configurations(nominal, perturbed, names)

params = {'Modulation_efficiency', 'Ratio'};
param_names = ["Modulation efficiency", "Nulling ratio"];
nParams = numel(params);
nConfig = numel(names);

path = "/Users/francesco/Library/CloudStorage/OneDrive-Personale/Università/SEMESTRE 12/TESI/Thesis/paperino/67f92ac5d9ea1e6dd82b036e/images/results/";

figure; hold on;

for i = 1:nConfig

    % Extract the simulation data for the i-th configuration
    simData = perturbed.(params{2})(i, :); 
    
    % Create a boxplot for this configuration; use a small offset for grouping
    boxplot(simData, 'positions', i, 'colors', 'b', 'symbol', 'r+');
    
    % Plot the nominal value as a horizontal line at the configuration's index
    nominalValue = nominal.(params{2})(i);
    plot([i-0.3, i+0.3], [nominalValue, nominalValue], 'r--','LineWidth',1.5);

end

set(gca, 'XTick', 1:nConfig, 'XTickLabel', names);
xlabel('Array Configuration');
ylabel(param_names{2});
title(['Distribution of ', params{2}]);
grid minor;
set(gca, "YScale", "log");
hold off;
    
export_figures("font_size", 10, "width", 20, "height", 12, "name", path + "_config_comparison_nulling_ratio")


figure; hold on;

violinplot(perturbed.(params{1})');

for i = 1:nConfig
    % Overlay the nominal value
    nominalValue = nominal.(params{1})(i);
    plot(i, nominalValue, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
end

set(gca, 'XTick', categorical(1:nConfig), 'XTickLabel', names);
xlabel('Array Configuration');

ylabel(param_names{1});
title(['Violin plot of ', param_names{1}]);
grid minor;
hold off;

export_figures("font_size", 10, "width", 20, "height", 12, "name", path + "_config_comparison_efficiency")


h = [figure(); figure()];
for i = 1:nConfig
    for p = 1:nParams

        figure(h(p));
        subplot(nConfig, 1, i)
        simData = perturbed.(params{p})(i, :);  % adjust indexing
        plot(simData, 'b.-');
        hold on
        nominalValue = nominal.(params{p})(i);
        yline(nominalValue, 'r--', "LineWidth", 1.5);
        xlabel('Simulation Index');
        ylabel(param_names{p});
        grid minor;
        if p == 1
            ylabel([names(i) ' - ' param_names{p}]);
        elseif p == 2
            set(gca, "YScale", "log")
        end
        hold off
    end
end

figure(h(1))
export_figures("font_size", 10, "width", 10, "height", 26, "name", path + "_config_comparison1")

figure(h(2))
export_figures("font_size", 10, "width", 10, "height", 26, "name", path + "_config_comparison2")

end