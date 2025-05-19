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
%   universes_to_plot[1]How many plots should be generated (default: 1)
%   plot_heatmap[bool]  If true (def.), gen. the heatmap or the histogram.
%   HZ_table[table]     Table that, for each stellar type, contains:
%       Stype[Nt x 1]   Type of the host star (like 'F', 'G', 'K', 'M').
%       HZ_min_10[Ntx1] Minimum HZ 10% percentile [AU]
%       HZ_min_90[Ntx1] Minimum HZ 90% percentile [AU]
%       HZ_min_mean[Ntx1] Minimum HZ mean value [AU]
%       ASHZ_min_10[Ntx1] Minimum HZ 10% percentile [rad]
%       ASHZ_min_90[Ntx1] Minimum HZ 90% percentile [rad]
%       ASHZ_min_mean[Ntx1] Minimum HZ mean value [rad]
%       HZ_max_10[Ntx1] Maximum HZ 10% percentile [AU]
%       HZ_max_90[Ntx1] Maximum HZ 90% percentile [AU]
%       HZ_max_mean[Ntx1] Maximum HZ mean value [AU]
%       ASHZ_max_10[Ntx1] Maximum HZ 10% percentile [rad]
%       ASHZ_max_90[Ntx1] Maximum HZ 90% percentile [rad]
%       ASHZ_max_mean[Ntx1] Maximum HZ mean value [rad]
%                       Default: empty.
%   filter_types[strings]Array of types that should be plotted.
%                       Default: ["F", "G", "K", "M"].
%   legend_position[string] Position of the legend (refer to legend docs).
%   show_legend[string] Either "show" or "off". 
%   embedded[struct]    All the setting for export, all in a single struct. 
%
% VERSION HISTORY:
%   2025-04-07 -------- 1.0
%   2025-04-08 -------- 1.1
%                     - Added several new plots to the analysis.
%   2025-04-10 -------- 1.2
%                     - Added option for universes_to_plot to reduce plots.
%   2025-04-24 -------- 1.2.1
%                     - Fixed exporting checks.
%   2025-05-08 -------- 1.3
%                     - Removed unused useless plots.
%                     - Small improvements to existing plots.
%                     - Below 3 simulations, use normal histogram.
%   2025-05-19 -------- 1.4
%                     - Introduced input parsing (and fixed it).
%                     - Use already defined rad2mas function.
%                     - New plots for HZ zones.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    exotable 
    IWA 
    OWA 
    ratio 
    export_setup.universes_to_plot = 1
    export_setup.embedded = struct()
    export_setup.plot_heatmap = true
    export_setup.HZ_table = table
    export_setup.legend_position = "best"
    export_setup.show_legend = "show"
    export_setup.filter_types = ["F", "G", "K", "M"]
end

% Compatibility mode
if ~isempty(fieldnames(export_setup.embedded))
    export_setup = export_setup.embedded;
end

% Filter out types
idx = ismember(exotable.Stype, export_setup.filter_types);
alltable = exotable; % Important for plotting purposes when filtered
exotable = exotable(idx, :);
idx = ismember(export_setup.HZ_table.Stype, export_setup.filter_types);
HZ_table = export_setup.HZ_table(idx, :);

universes = unique(exotable.Nuniverse);
Nu = length(universes);
Ns = size(exotable.yields, 2);
markers = {'o','s','^','d','v','>','<','p','h'};
unique_types = unique(exotable.Stype);
unique_alltypes = unique(alltable.Stype);

cols = get_colours(2);

style_colors;

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
        idx = find(ismember(unique_alltypes, current_type) == 1);
        
        if strcmp(current_type, "")
            current_type = "N.A.";
        end
        
        if length(unique_types) < length(markers)
            % Plot all points in blue
            scatter(rad2mas(T.ang_sep_rad(type_rows)), T.contrast(type_rows), ...
                    [], cols(2, :), markers{idx}, 'filled', 'DisplayName', current_type);
        else
            scatter(rad2mas(T.ang_sep_rad(type_rows)), T.contrast(type_rows), ...
                    [], 'filled', 'DisplayName', current_type);
        end
    end

    for j = 1:length(unique_types)
        current_type = unique_types{j};
        type_rows = strcmp(T.Stype, current_type);
        idx = find(ismember(unique_alltypes, current_type) == 1);
        
        if strcmp(current_type, "")
            current_type = 'N.A.';
        end
        
        if length(unique_types) < length(markers)
            % Overlay resulting points
            yield_one_rows = max(type_rows & (T.yields == 1), [], 2);
            scatter(rad2mas(T.ang_sep_rad(yield_one_rows)), T.contrast(yield_one_rows), ...
                    [], cols(1, :), markers{idx}, 'filled', 'DisplayName', [current_type ' (Detected)']);
        else
            % Overlay resulting points
            yield_one_rows = max(type_rows & (T.yields == 1), [], 2);
            scatter(rad2mas(T.ang_sep_rad(yield_one_rows)), T.contrast(yield_one_rows), ...
                    [], cols(1, :), 'filled', 'DisplayName', [current_type ' (Detected)']);
        end
    end

    xline(rad2mas(IWA), '--', "DisplayName", "IWA", "Color", ui_colours(1, :), "LineWidth", 1.5)
    xline(rad2mas(OWA), '--', "DisplayName", "OWA", "Color", ui_colours(2, :), "LineWidth", 1.5)
    yline(ratio, '--', "DisplayName", "Best nulling ratio", "Color", ui_colours(end, :), "LineWidth", 1.5)

    set(gca, 'XScale', 'log', 'YScale', 'log'); 
    legend(export_setup.show_legend, 'Location', export_setup.legend_position, 'NumColumns', 3);
    xlabel("Angular separation [mas]")
    ylabel("Flux contrast at observation wavelength")
    grid minor;
    
    if ~isempty(HZ_table)

        % Get the current y-axis limits to restore them later
        ylims = get(gca, 'YLim');
        
        % Only for HZ plots, extract some colours
        Hcols = get_colours(size(export_setup.HZ_table, 1), "ui");

        for j = 1:size(HZ_table, 1)
            
            % Finish colour extraction (same colour for type)
            idx = find(ismember(export_setup.HZ_table.Stype, export_setup.filter_types(j)) == 1);

            x_min_mean = rad2mas(HZ_table.ASHZ_min_mean(j));
            x_max_mean = rad2mas(HZ_table.ASHZ_max_mean(j));
            x_min_10 = rad2mas(HZ_table.ASHZ_min_10(j));
            x_max_90 = rad2mas(HZ_table.ASHZ_max_90(j));
        
            % Plot the first rectangle (alpha 0.5, no border)
            rectangle('Position', [x_min_mean, 1e-15, x_max_mean - x_min_mean, 1e10], ...
                'FaceColor', Hcols(idx, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5);

            % Faking line for legends (...)
            plot(NaN,NaN,'LineWidth',10,'LineStyle','-','Color',Hcols(idx,:), ...
                "DisplayName", sprintf("Mean HZ for type %s", export_setup.filter_types(j)));
            
            % Plot the second rectangle (no background color, thick border)
            rectangle('Position', [x_min_10, 1e-15, x_max_90 - x_min_10, 1e10], ...
                'FaceColor', 'none', 'EdgeColor', Hcols(idx, :), 'LineWidth', 2, 'LineStyle', ":");
           
            % Faking line for legends (...)
            plot(NaN,NaN,'LineWidth',1.5,'LineStyle',':','Color',Hcols(idx,:), ...
                "DisplayName", sprintf("External percentiles HZ for type %s", export_setup.filter_types(j)));

        end

        % Reapply the original y-axis limits
        set(gca, 'YLim', ylims);
    end

    if isstruct(export_setup) && isfield(export_setup, "figures") && isstruct(export_setup.figures.yields)
        export_figures("embedded", export_setup.figures.yields)
    end
end


if export_setup.plot_heatmap

    yields = exotable.yields;           % N_planets x Ns (1: detected, 0: no)       
    universes = exotable.Nuniverse;     % N_planets x 1  (0..0, 1..1, ...)  
    
    total_yield_matrix = zeros(Nu, Ns); % Nu x Ns
    
    for u = 1:Nu
        % Logical indexing for current universe and sum to get plot
        idx = universes == u-1;              
        total_yield_matrix(u, :) = sum(yields(idx, :));
    end
    
    if Ns > 3
    
        figure;
        imagesc(total_yield_matrix);
        xlabel('Simulation');
        ylabel('Universe');
        title('Total Yields');
        grid minor;
        colorbar;
        clim([0, max(total_yield_matrix, [], "all")])
        colormap(darkBlue);
    
    else
        
        figure;
        bar(total_yield_matrix, 'grouped'); 
        xlabel('Universe');
        ylabel('Total Yields');
        if Ns > 1
            legend(arrayfun(@(x) sprintf('Simulation %d', x), 1:Ns, 'UniformOutput', false));
        end
        grid on;
    
    end
    
    if isstruct(export_setup) && isfield(export_setup, "figures") && isstruct(export_setup.figures.total_yield)
        export_figures("embedded", export_setup.figures.total_yield)
    end
end


end