function export_sensitivity_table(tableData, names, filename)
%EXPORT_SENSITIVITY_TABLE Export a LaTeX tabular with the results of the
%sensitivity analysis to the baseline and the aperture ratios, as found
%from plot_array_size_sensitivity.
%
% INPUTS:
%   tableData[table]    Table or struct with fields
%       - B_opt[Nx3]    Optimal baselines for weights 0.4, 0.5, 0.6 [m]
%       - AR_opt[Nx3]   Optimal aperture ratio for weights 0.4, 0.5, 0.6
%       - NR_opt[Nx6]   Nulling ratios columns:
%                           1-2 (p=0.4), 3-4 (p=0.5), 5-6 (p=0.6)
%       - IWA_opt[Nx6]  IWAs in rad (as for NR)
%       - detections[Nx6] Detected planets (as for NR)
%   names[strings]      Missions names
%   filename[string]    Filename (without .tex)
%
% VERSION HISTORY:
%   2025-05-15 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename + ".tex", 'w');


% Header
fprintf(fid, '    \\begin{tabular}{lccccccccc}\n');
fprintf(fid, '        \\toprule\n');
fprintf(fid, '        & & \\multicolumn{2}{c}{\\textbf{Optimisation parameter}} & \\multicolumn{2}{c}{\\textbf{Nulling ratio}} & \\multicolumn{2}{c}{\\textbf{Inner working angles [mas]}} & \\multicolumn{2}{c}{\\textbf{Detections}} \\\\\n');
fprintf(fid, '        \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}\n');
fprintf(fid, '        \\textbf{Mission} & \\textbf{Weights} & Baseline [m] & Aperture ratio & Baseline [m] & Aperture ratio & Baseline [m] & Aperture ratio & Baseline [m] & Aperture ratio \\\\\n');
fprintf(fid, '        \\midrule \\midrule\n');

% Loop
for i = 1:length(names)
    for w = 1:3
        if w == 1
            fprintf(fid, '        \\multirow{3}{*}{%s} & 0.4', names{i});
        elseif w == 2
            fprintf(fid, '         & 0.5');
        else
            fprintf(fid, '         & 0.6');
        end

        % Optimisation
        b = tableData.B_opt(i, w);
        ar = tableData.AR_opt(i, w);

        % Nulling ratio
        nr_b = log10(tableData.NR_opt(i, (w-1)*2 + 1));
        nr_ar = log10(tableData.NR_opt(i, (w-1)*2 + 2));

        % IWA
        iwa_b = rad2mas(tableData.IWA_opt(i, (w-1)*2 + 1));
        iwa_ar = rad2mas(tableData.IWA_opt(i, (w-1)*2 + 2));

        % Detections
        det_b = tableData.detections(i, (w-1)*2 + 1);
        det_ar = tableData.detections(i, (w-1)*2 + 2);

        % Write
        fprintf(fid, ' & %.2f & %.2f & $10^{%.2f}$ & $10^{%.2f}$ & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
            b, ar, nr_b, nr_ar, iwa_b, iwa_ar, det_b, det_ar);
    end

    % Empty row
    fprintf(fid, '        \\\\ \\midrule\n');
end

% Footer
fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');

fclose(fid);
end
