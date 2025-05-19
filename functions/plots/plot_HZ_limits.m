function plot_HZ_limits(starTable, summaryTable)
%PLOT_HZ_LIMITS Visualizes HZ bounds for different star types.
%
% OPTIONAL INPUTS:
%   startable[table]    Unique stars of the catalogue. With fields:
%       Nuniverse[Nsx1] Universe associated to that planet.
%       Nstar[Ns x 1]   Univoque number of the star in the universe.
%       Rs[Ns x 1]      Radius of the star. [R_Sun]
%       Ts[Ns x 1]      Effective temperature of star. [K]
%       Ds[Ns x 1]      Distance to Sun of host star. [pc]
%       Ms[Ns x 1]      Mass of the star. [M_Sun]
%       Stype[Ns x 1]   Type of the host star (like 'F', 'G', 'K', 'M').
%       RA[Ns x 1]      Right ascension of the star [deg]
%       Dec[Ns x 1]     Declination of the star [deg]
%       HZ_min[Ns x 1]  Minimum position of habitable zone [AU]
%       HZ_max[Ns x 1]  Maximum position of habitable zone [AU]
%       AngSep_HZ_min[Nsx1] Min pos of HZ in ang. sep. from Earth [rad]
%       AngSep_HZ_max[Nsx1] Max pos of HZ in ang. sep. from Earth [rad]
%
%   summaryTable[table] Table that, for each stellar type, contains:
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
%
% NOTES:
%   - The function can be called with just starTable, just summaryTable or
%     both inputs, in order to visualise different degrees of results. 
%
% REFERENCES:
%   Kopparapu RK, Ramirez R, Kasting JF, Eymet V, Robinson TD, Mahadevan S,
%   et al. Habitable zones around main-sequence stars: new estimates. ApJ.
%   2013 Feb 26;765(2):131.
%
% VERSION HISTORY:
%   2025-05-19 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
grid minor;
box on;
style_colors;

% Unique spectral types and color palette
if istable(starTable)
    stypes = unique(starTable.Stype);
    cols = get_colours(length(stypes), "ui");
end

if nargin > 1 && ~isempty(summaryTable)
  
    % Get unique spectral types from summary table
    summaryStypes = unique(summaryTable.Stype);
    cols = get_colours(length(summaryStypes), "ui");
    
    % Plot rectangles for each spectral type in summary table
    for i = 1:numel(summaryStypes)
        thisType = summaryStypes(i);
        idx = strcmp(summaryTable.Stype, thisType);
        
        % Get values
        x_min_mean = summaryTable.HZ_min_mean(idx);
        x_max_mean = summaryTable.HZ_max_mean(idx);
        x_min_10 = summaryTable.HZ_min_10(idx);
        x_max_90 = summaryTable.HZ_max_90(idx);
        
        y_min_mean = summaryTable.ASHZ_min_mean(idx);
        y_max_mean = summaryTable.ASHZ_max_mean(idx);
        y_min_10 = summaryTable.ASHZ_min_10(idx);
        y_max_90 = summaryTable.ASHZ_max_90(idx);
        
        % Plot inner rectangle with solid color
        rectangle('Position', [x_min_mean, rad2mas(y_min_mean), x_max_mean - x_min_mean, rad2mas(y_max_mean - y_min_mean)], ...
            'FaceColor', cols(i,:), 'EdgeColor', cols(i,:));
        
        % Plot outer rectangle with 0.5 alpha
        rectangle('Position', [x_min_10, rad2mas(y_min_10), x_max_90 - x_min_10, rad2mas(y_max_90 - y_min_10)], ...
            'FaceColor', cols(i,:), 'EdgeColor', cols(i,:), 'FaceAlpha', 0.5);
        
        % Faking line for legends (...)
        line(NaN,NaN,'LineWidth',10,'LineStyle','-','Color',cols(i,:));

        if isempty(starTable)
            legendEntries{i} = sprintf('Spectral type %s', cell2mat(thisType));
        end

    end
end

if istable(starTable)
    legendEntries = cell(numel(stypes),1);
    
    % Plot ovals per spectral type
    for i = 1:numel(stypes)
        thisType = stypes(i);
        idx = strcmp(starTable.Stype, thisType);
    
        % Get values
        x = (starTable.HZ_min(idx) + starTable.HZ_max(idx)) / 2;
        dx = (starTable.HZ_max(idx) - starTable.HZ_min(idx)) / 2;
    
        y = (starTable.AngSep_HZ_min(idx) + starTable.AngSep_HZ_max(idx)) / 2;
        dy = (starTable.AngSep_HZ_max(idx) - starTable.AngSep_HZ_min(idx)) / 2;
    
        % Plot each "oval" as an error bar
        errorbar(x, rad2mas(y), rad2mas(dy), rad2mas(dy), dx, dx, ...
            'o', 'MarkerSize', 5, ...
            'MarkerEdgeColor', cols(i,:), ...
            'MarkerFaceColor', cols(i,:), ...
            'Color', cols(i,:), ...
            'LineWidth', 1, ...
            'LineStyle', 'none');
    
        legendEntries{i} = sprintf('Spectral type %s', cell2mat(thisType));
    end
end

xlabel('Habitable Zone distance [AU]');
ylabel('Angular separation of HZ [mas]');
title('Habitable Zones for Different Spectral Types');
legend(legendEntries, 'Location', 'best');
set(gca, "YScale", "log")
hold off;
end
