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
%   2025-04-24 -------- 1.2.1
%                     - Fixed exporting checks.
%   2025-05-08 -------- 1.3
%                     - Removed unused useless plots.
%                     - Small improvements to existing plots.
%                     - Below 3 simulations, use normal histogram.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

universes = unique(exotable.Nuniverse);
Nu = length(universes);
Ns = size(exotable.yields, 2);
rad2mas = ((pi/180)/3600 / 1e3)^(-1);     
markers = {'o','s','^','d','v','>','<','p','h'};
unique_types = unique(exotable.Stype);

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
        
        if strcmp(current_type, "")
            current_type = "N.A.";
        end
        
        if length(unique_types) < length(markers)
            % Plot all points in blue
            scatter(T.ang_sep_rad(type_rows) * rad2mas, T.contrast(type_rows), ...
                    [], cols(2, :), markers{j}, 'filled', 'DisplayName', current_type);
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
                    [], cols(1, :), markers{j}, 'filled', 'DisplayName', [current_type ' (Detected)']);
        else
            % Overlay resulting points
            yield_one_rows = max(type_rows & (T.yields == 1), [], 2);
            scatter(T.ang_sep_rad(yield_one_rows) * rad2mas, T.contrast(yield_one_rows), ...
                    [], cols(1, :), 'filled', 'DisplayName', [current_type ' (Detected)']);
        end
    end

    xline(IWA * rad2mas, '--', "DisplayName", "IWA", "Color", ui_colours(1, :), "LineWidth", 1.5)
    xline(OWA * rad2mas, '--', "DisplayName", "OWA", "Color", ui_colours(2, :), "LineWidth", 1.5)
    yline(ratio, '--', "DisplayName", "Best nulling ratio", "Color", ui_colours(end, :), "LineWidth", 1.5)

    set(gca, 'XScale', 'log', 'YScale', 'log'); 
    legend('show', 'Location', 'best', 'NumColumns', 3);
    xlabel("Angular separation [mas]")
    ylabel("Flux contrast at observation wavelength")
    grid minor;

    if isstruct(export_setup) && isfield(export_setup, "figures") && isstruct(export_setup.figures.yields)
        export_figures("embedded", export_setup.figures.yields)
    end
end

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