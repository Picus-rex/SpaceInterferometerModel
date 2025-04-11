function plot_ppop_yield(exotable, IWA, OWA, ratio, export_setup)
%PLOT_PPOP_YIELD Create plots from the results of get_ppop_yield
%
% INPUTS:
%   exotable[table]     Resulting table from get_ppop_yield, filled.
%   IWA[1]              (Best) inner working angle. [rad]
%   OWA[1]              (Best) outer working angle. [rad]
%   ratio[1]            (Best) nulling ratio of the system. [-]
%
% ARGUMENT OPTIONAL INPUTS:
%   export_setting[struct] Setting for export, all in a single struct. 
%
% VERSION HISTORY:
%   2025-04-07 -------- 1.0
%   2025-04-08 -------- 1.1
%                     - Added several new plots to the analysis.
%   2025-04-10 -------- 1.2
%                     - Added option for universes_to_plot to reduce plots.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

universes = unique(exotable.Nuniverse);
Nu = length(universes);
Ns = size(exotable.yields, 2);
rad2mas = ((pi/180)/3600 / 1e3)^(-1);     
markers = {'o','s','^','d','v','>','<','p','h'};
unique_types = unique(exotable.Stype);

if isfield(export_setup, "universes_to_plot")
    selection = randi(Nu, export_setup.universes_to_plot);
else
    selection = randi(Nu, 5);
end
selection = selection(1, :);

for i = 1:length(selection)

    T = exotable(exotable.Nuniverse == universes(selection(i)), :);

    figure; hold on;

    % Loop through each unique type and plot with the corresponding marker
    for j = 1:length(unique_types)
        current_type = unique_types{j};
        type_rows = strcmp(T.Stype, current_type);
        
        if strcmp(current_type, "")
            current_type = "N.A.";
        end
        
        if length(unique_types) < length(markers)
            % Plot all points in blue
            scatter(T.ang_sep_rad(type_rows) * rad2mas, T.contrast(type_rows), ...
                    [], 'b', markers{j}, 'filled', 'DisplayName', current_type);
        else
            scatter(T.ang_sep_rad(type_rows) * rad2mas, T.contrast(type_rows), ...
                    [], 'filled', 'DisplayName', current_type);
        end
    end

    for j = 1:length(unique_types)
        current_type = unique_types{j};
        type_rows = strcmp(T.Stype, current_type);
        
        if strcmp(current_type, "")
            current_type = 'N.A.';
        end
        
        if length(unique_types) < length(markers)
            % Overlay resulting points
            yield_one_rows = max(type_rows & (T.yields == 1), [], 2);
            scatter(T.ang_sep_rad(yield_one_rows) * rad2mas, T.contrast(yield_one_rows), ...
                    [], 'r', markers{j}, 'filled', 'DisplayName', [current_type ' (Detected)']);
        else
            % Overlay resulting points
            yield_one_rows = max(type_rows & (T.yields == 1), [], 2);
            scatter(T.ang_sep_rad(yield_one_rows) * rad2mas, T.contrast(yield_one_rows), ...
                    [], 'r', 'filled', 'DisplayName', [current_type ' (Detected)']);
        end
    end

    xline(IWA * rad2mas, '--', "DisplayName", "IWA (Mean value)")
    xline(OWA * rad2mas, '--', "DisplayName", "OWA")
    yline(ratio, '--', "DisplayName", "Best achievable nulling ratio")

    set(gca, 'XScale', 'log', 'YScale', 'log'); 
    legend('show', 'Location', 'best', 'NumColumns', 3);
    xlabel("Angular separation [mas]")
    ylabel("Flux contrast at observation wavelength")
    grid minor;

    if isstruct(export_setup) && isfield(export_setup, "figures") && export_setup.figures.yields
        export_figures("embedded", export_settings.figures.yields)
    end
end

yields = exotable.yields;        
universes = exotable.Nuniverse;  

mean_yield_matrix = zeros(Nu, Ns); % Nu x Ns
total_yield_matrix = zeros(Nu, Ns); % Nu x Ns

for u = 1:Nu
    idx = universes == u-1;              % Logical indexing for current universe
    mean_yield_matrix(u, :) = mean(yields(idx, :), 1); % Mean across planets in universe
    total_yield_matrix(u, :) = sum(yields(idx, :));
end

figure;
imagesc(mean_yield_matrix);
xlabel('Simulation');
ylabel('Universe');
title('Mean Yields');
grid minor;
colorbar;
if isstruct(export_setup) && isfield(export_setup, "figures") && export_setup.figures.mean_yield
    export_figures("embedded", export_settings.figures.mean_yield)
end

figure;
imagesc(total_yield_matrix);
xlabel('Simulation');
ylabel('Universe');
title('Total Yields');
grid minor;
colorbar;
if isstruct(export_setup) && isfield(export_setup, "figures") && export_setup.figures.total_yield
    export_figures("embedded", export_settings.figures.total_yield)
end

if Nu > 1
    figure;
    boxplot(mean_yield_matrix, 'Labels', string(1:Ns));
    xlabel('Simulation Index');
    ylabel('Mean Yield Across Universes');
    title('Boxplot of Mean Yields per Simulation');
    grid minor;
    if isstruct(export_setup) && isfield(export_setup, "figures") && export_setup.figures.boxplot
        export_figures("embedded", export_settings.figures.boxplot)
    end
end

figure;
hold on;

for u = 1:Nu
    plot(mean_yield_matrix(u, :), '-o', 'DisplayName', ['Universe ' num2str(u)]);
end

xlabel('Simulation Index');
ylabel('Mean Yield');
title('Yield Trend Across Simulations per Universe');
legend('Location', 'best');
grid minor;
if isstruct(export_setup) && isfield(export_setup, "figures") && export_setup.figures.trend
    export_figures("embedded", export_settings.figures.trend)
end

end