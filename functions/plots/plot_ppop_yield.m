function plot_ppop_yield(exotable, IWA, OWA, export_setup)

universes = unique(exotable.Nuniverse);
Nu = length(universes);
arcsec2rad = (pi/180)/3600;     % arcsec to radians conversion factor
markers = {'o','s','^','d','v','>','<','p','h'};
unique_types = unique(exotable.Stype);

for i = 1:Nu

    T = exotable(exotable.Nuniverse == universes(i), :);

    figure; hold on;

    % Loop through each unique type and plot with the corresponding marker
    for j = 1:length(unique_types)
        current_type = unique_types{j};
        type_rows = strcmp(T.Stype, current_type);
        
        % Plot all points in blue
        scatter(T.ang_sep_rad(type_rows) / arcsec2rad, T.contrast(type_rows), ...
                [], 'b', markers{j}, 'filled', 'DisplayName', current_type);
    end

    for j = 1:length(unique_types)
        current_type = unique_types{j};
        type_rows = strcmp(T.Stype, current_type);

        % Overlay resulting points
        yield_one_rows = max(type_rows & (T.yields == 1), [], 2);
        scatter(T.ang_sep_rad(yield_one_rows) / arcsec2rad, T.contrast(yield_one_rows), ...
                [], 'r', markers{j}, 'filled', 'DisplayName', [current_type ' (Detected)']);
    end

    xline(IWA / arcsec2rad, '--', "DisplayName", "IWA (Mean value)")
    xline(OWA / arcsec2rad, '--', "DisplayName", "OWA")

    set(gca, 'XScale', 'log', 'YScale', 'log'); 
    legend('show', 'Location', 'best', 'NumColumns', 2);
    xlabel("Angular separation [mas]")
    ylabel("Flux contrast at observation wavelength")
    grid minor;
end



end