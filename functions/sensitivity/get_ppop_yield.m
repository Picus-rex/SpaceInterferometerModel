function [yields, candidates] = get_ppop_yield(IWAs, OWA, ratios, lambda)

Ns = size(IWAs, 1);

opts = detectImportOptions('others/TestPlanetPopulation.txt', 'NumHeaderLines', 1);
exotable = readtable('others/TestPlanetPopulation.txt', opts);

% Flux contrast
AU = 1.496e11;                  % Astronomical unit in meters
h = 6.626e-34;                  % Planck constant [J s]
c = 3e8;                        % Speed of light [m/s]
k = 1.381e-23;                  % Boltzmann constant [J/K]
R_earth = 6.371e6;              % Earth radius [m]
R_sun   = 6.957e8;              % Solar radius [m]

for i = 1:10

    universeToPlot = i-1;
    rows = exotable.Nuniverse == universeToPlot;
    T = exotable(rows, :);
    
    
    arcsec2rad = (pi/180)/3600;
    ang_sep_rad(1:length(T.AngSep), i) = T.AngSep * arcsec2rad;
    
    
    B = @(T_val) (2*h*c^2 ./ lambda^5) ./ (exp(h*c./(lambda*k*T_val)) - 1);

    thermal_contrast = ((T.rp * R_earth).^2 .* B(T.Tp)) ./ ((T.Rs * R_sun).^2 .* B(T.Ts));
    reflected_contrast = T.AgeomMIR .* ((T.rp * R_earth) ./ (T.ap * AU)).^2;

    contrast(1:length(T.AngSep), i) = thermal_contrast + reflected_contrast;

end

Nu = size(contrast, 2);
yields = zeros(Nu, Ns);
candidates = zeros(Nu, Ns);

for i = 1:Ns
    
    IWA = IWAs(i);
    ratio = rms(ratios(:, i));

    for j = 1:Nu
        for k = 1:size(contrast, 1)
            
            if ang_sep_rad(k, j) > IWA && ang_sep_rad(k, j) < OWA
                
                candidates(j, i) = candidates(j, i) + 1;

            end

            if ang_sep_rad(k, j) > IWA && ang_sep_rad(k, j) < OWA && contrast(k, j) > ratio % * 1e-6
                
                yields(j, i) = yields(j, i) + 1;
                
            end
        end
    end

end



end








% figure;
% hold on; 
% set(gca, 'XScale', 'log', 'YScale', 'log'); 
% 
% starTypes = unique(T.Stype);
% 
% % Define markers and colors for different groups (adjust as desired)
% markers = {'o','s','^','d','v','>','<','p','h'};
% %colors = lines(numel(starTypes));
% 
% style_colors;
% 
% for i = 1:length(starTypes)
% 
%     % Select indices for the current star type
%     idx = strcmp(T.Stype(:), {starTypes{i}});
% 
%     % Create a log–log scatter plot for this group
%     scatter(ang_sep_rad(idx) * 648000000/pi, contrast(idx), 40, ...
%         colours(i,:), 'filled', 'DisplayName', starTypes{i});
% 
% end
% 
% % xline(resol, '--', "LineWidth", 1.5, "DisplayName", "Angular resolution");
% % xline(IWA, 'b--', "LineWidth", 1.5, "DisplayName", "IWA");
% % xline(OWA, 'r--', "LineWidth", 1.5, "DisplayName", "OWA");
% 
% xlabel('Angular Separation [mas]');
% ylabel('Flux Contrast at 10 \mum');
% % title('Exoplanet Population: Angular Separation vs. Flux Contrast');
% legend('show','Location','best');
% grid on;
% hold off;
% 
% 
% 
% 
% 
% 
% end