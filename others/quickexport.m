% Loop through all open figure handles
figs = findall(0, 'Type', 'figure');

for i = 1:length(figs)
    figure(figs(i)); % Make the figure current

    % Set font for all objects in the figure
    set(findall(figs(i), '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(findall(figs(i), '-property', 'FontSize'), 'FontSize', 10);

    % Remove title from current axes
    ax = gca;
    title(ax, "");

    % Set figure and paper size
    set(figs(i), 'Units', 'centimeters', 'Position', [0, 0, 10, 5]);
    set(figs(i), 'PaperUnits', 'centimeters', 'PaperSize', [10, 5]);
    set(figs(i), 'PaperPositionMode', 'auto');

    % Save as PDF with high resolution
    print(figs(i), sprintf("exports/images/validation/%d", i), '-dpdf', '-r300');
end