function plot_variation_map(Ns, theta_x, theta_y, Maps, Nominal_Map, export_setup)
%PLOT_VARIATION_MAP Generate a series of plots to study how variations maps
%changes due to perturbations. This function is only meant to be called
%from perturbate_system. 

% Conversion factor and colours where needed
style_colors;
rad2mas = 1e3 * (3600 * 180) / pi;

% FIGURES 2
% Extract random maps
permute_data = randperm(Ns);
for i = 1:export_setup.perturbed_map_plotting_number
    plot_transmission_map(theta_x(1, :), Maps(:, :, permute_data(i)));
end

% FIGURE 3
% Compute the difference maps elements-wise without doing a for, then
% calculate the standard deviation (variability) across simulations
diffMaps = bsxfun(@minus, Maps, Nominal_Map);
stdMap = std(diffMaps, 0, 3);  % 0 indicates normalization by N-1

figure;
%contourf(theta_x, theta_y, stdMap, 20, 'LineColor', 'none'); 
contourf(theta_x*rad2mas, theta_y*rad2mas, stdMap, 20); 
colormap(darkBlue); colorbar;
title('Standard Deviation Map of Transmission Variability');
xlabel('\theta_x [mas]');
ylabel('\theta_y [mas]');
axis square;

% Determine the size of the maps
[m, n, ~] = size(Maps);

if Ns > 1
    % FIGURE 4
    % Each simulation is vectorized to form a row vector.
    dataMatrix = reshape(Maps, [m*n, Ns])';  % dimensions: Ns x (m*n)
    
    % Remove the mean from each pixel across simulations and erform PCA on 
    % the centered data, then reshape the first principal mode into 2D map
    meanMapVec = mean(dataMatrix, 1);
    dataMatrix_centered = dataMatrix - meanMapVec;
    [coeff, score, ~] = pca(dataMatrix_centered);
    pc1 = reshape(coeff(:,1), [m, n]);
    
    % Plot the first principal component mode as a contour map
    figure;
    contourf(theta_x*rad2mas, theta_y*rad2mas, pc1, 20);
    colormap(darkBlue); colorbar;
    title('First Principal Component Mode');
    xlabel('\theta_x [mas]');
    ylabel('\theta_y [mas]');
    axis square;
    
    % FIGURE 5
    % Plot the PC1 scores for each simulation
    figure;
    plot(score(:,1), 'o-', "LineWidth", 1.5, "Color", colours(1, :));
    xlabel('Simulation Number');
    ylabel('PC1 Score');
    title('PC1 Scores Across Simulations');
    grid minor;
end

% FIGURE 6
% Plot a few representative difference maps
figure;
nPlots = min(4, Ns);  % choose up to 4 maps to display
for i = 1:nPlots
    subplot(2,2,i);
    contourf(theta_x*rad2mas, theta_y*rad2mas, diffMaps(:,:,permute_data(i)), 20);
    colormap(darkBlue); colorbar;
    title(['Difference Map, Simulation ', num2str(permute_data(i))]);
    xlabel('\theta_x [mas]');
    ylabel('\theta_y [mas]');
    axis square;
end

% FIGURE 7
% Plot the cumulative distribution for a selected pixel
% Choose pixel coordinates, e.g., the center of the map
ix = round(m/2);
iy = round(n/2);
pixelDiffs = squeeze(diffMaps(ix, iy, :));  % differences at the selected pixel

% Sort the differences and calculate the empirical CDF
[sortedDiffs, ~] = sort(pixelDiffs);
cdfValues = (1:Ns) / Ns;

figure;
plot(sortedDiffs, cdfValues, 'o-', "LineWidth", 1.5, "Color", colours(1, :));
xlabel('Difference Value at Pixel');
ylabel('Cumulative Distribution');
title(['CDF of Differences at Pixel (', num2str(ix), ', ', num2str(iy), ')']);
grid minor;

end