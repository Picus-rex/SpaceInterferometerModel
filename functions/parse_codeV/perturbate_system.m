function [ratios, Maps] = perturbate_system(N, amplitudes_nominal, ...
    phases_nominal, positions_nominal, theta_star, lambda, ...
    phases_perturbed, combinations, theta_x, theta_y, export_setup)
%PERTURBATE_SYSTEM From the nominal characteristics of the system, computes
%the variations based on the provided phases shift. 
%
% INPUTS:
%   N[1]                Number of apertures [-]
%   amplitude_nom[Nx1]  Vector of nominal aperture amplitudes [m^2]
%   phases_nom[Nx1]     Vector of nominal phase shifts for each aperture
%   positions_nom[Nx2]  Matrix of (x, y) nominal positions of the apertures
%   theta_star[1]       Angular dimension of the star [rad]
%   lambda[1]           Wavelength [m]
%   phases_pert[NxNs]   Resulting perturbed phases for one aperture over
%                       all the rays (N points) over Ns simulations
%   combinations[Nx1]   Vector of beam combiner coefficients
%   theta_x, theta_y    Meshgrid of angular coordinates
% 
% ARGUMENT INPUTS:
%   create_plots[bool]  If true, creates several plots.
%   type[string]        String following "%unit_%de(+-)%d",
%                       where the exponential number is reversed and
%                       recognised as a scale (if m_1e-6 is given, for
%                       example, the string is converted to micro m). 
%   perturbed_map_plotting_number[1] Number of perturbed maps to plot (must
%                       be below Ns.
%
% OUTPUTS:
%   ratios[N x Ns]      Nulling ratios associated to every provided OPD.
%   Maps[TX x TY x Ns]  3D matrix where the first two dimensions are the
%                       standard normalised response function along the
%                       angular coordinates provided.
%
% VERSION HISTORY:
%   2025-04-02 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    N {mustBeInteger}
    amplitudes_nominal (:, :) {mustBeNumeric}
    phases_nominal (:, :) {mustBeNumeric}
    positions_nominal (:, :) {mustBeNumeric}
    theta_star (1, 1) {mustBeNumeric}
    lambda (1, 1) {mustBeNumeric}
    phases_perturbed (:, :) {mustBeNumeric}
    combinations (:, :) {mustBeNumeric}
    theta_x = load_grid("x")
    theta_y = load_grid("y")
    export_setup.create_plots = true
    export_setup.type = "m_1e-6"
    export_setup.perturbed_map_plotting_number = 3
end

% Number of series, obtain OPDs as well as the nominal map
Ns = size(phases_perturbed, 2);
OPDs = phase2opd(phases_perturbed, lambda);
nom_table = compute_response_function("lambda", lambda, "N", N, ...
        "positions", positions_nominal, "A", amplitudes_nominal, ...
        "phase_shifts", phases_nominal, "combinations", combinations, ...
        "theta_x", theta_x, "theta_y", theta_y);
Nominal_Map = nom_table(:, :).T_standard;

% Allocation and reordering of the vector
phases_perturbed_all = reshape(phases_perturbed, [], 1);
ratios = zeros(length(phases_perturbed_all), 1);
perturbed_vect = zeros(size(phases_nominal));
Maps = zeros(length(theta_x), length(theta_y), Ns);

% Add the single perturbation to the phases to compute the nulling ratio
for i = 1:length(phases_perturbed_all)
    perturbed_vect(1) = phases_perturbed_all(i);
    phases = phases_nominal + perturbed_vect;

    amplitudes_modified = amplitudes_nominal .* combinations;

    ratios(i) = compute_nulling_ratio(N, amplitudes_modified, phases, ...
                                    positions_nominal, theta_star, lambda);
end

for i = 1:Ns
    perturbed_vect(1) = rms(phases_perturbed(:, i));
    phases = phases_nominal + perturbed_vect;

    map_table = compute_response_function("lambda", lambda, "N", N, ...
        "positions", positions_nominal, "A", amplitudes_nominal, ...
        "phase_shifts", phases, "combinations", combinations, "theta_x", ...
        theta_x, "theta_y", theta_y);

    Maps(:, :, i) = map_table(:, :).T_standard;

end

% Reformat the vector 
ratios = reshape(ratios, [], Ns);

if export_setup.create_plots
    
    % Conversion factor and colours where needed
    style_colors;
    rad2mas = 1e3 * (3600 * 180) / pi;

    % Get visualisation style
    [scale, scale_tag] = get_scale_plots(export_setup.type);
    elem_label = sprintf("OPD [%s]", scale_tag);
    
    % FIGURE 1
    figure; hold on;
    cols = get_colours(Ns);
    for i = 1:Ns
        scatter(OPDs(:, i) * scale, ratios(:, i), 36, cols(i, :), ...
                           "filled", "DisplayName", sprintf("Sim. %d", i));
    end
    xlabel(elem_label)
    ylabel("Nulling ratios")
    legend("Location","best");
    grid minor;
    set(gca, "YScale", "log");
    
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

    % FIGURE 4
    % Determine the size of the maps
    [m, n, ~] = size(Maps);
    
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

    % FIGURE 8
    maps_named = table();
    names_plot = strings(Ns, 1);
    for i = 1:Ns
        names_plot(i) = sprintf("Simulation %d", i);
        maps_named.(names_plot(i)) = Maps(:, :, i);
    end
    plot_transmission_map_monodirectional(theta_x(1, :), maps_named, names_plot);

end



end