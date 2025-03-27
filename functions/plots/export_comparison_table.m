function export_comparison_table(tableData, names, filename)

% Open file for writing
fid = fopen(filename + ".tex", 'w');

% Write the LaTeX table preamble
fprintf(fid, '\\begin{table}[htb!]\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\small\n');
fprintf(fid, '    \\begin{tabular}{lcccccccc}\n');
fprintf(fid, '        \\toprule\n');
fprintf(fid, '         & &  \\multicolumn{2}{c}{\\textbf{Position}} \\\\ \\cmidrule(lr){3-4}\n');
fprintf(fid, '        \\textbf{Configuration} & \\textbf{Ap.} & \\textbf{X}$_j$ [L] & \\textbf{Y}$_j$ [L] & $\\boldsymbol{\\varphi}$ [deg] & $\\boldsymbol{\\eta}_{\\textbf{mod,}\\boldsymbol{\\infty}}$ & $\\boldsymbol{\\vartheta}_\\textbf{IWA}$ $\\left[ \\nicefrac{\\lambda}{L} \\right]$ & $\\boldsymbol{\\vartheta}_\\textbf{FWHM}$ $\\left[ \\nicefrac{\\lambda}{L} \\right]$ & \\textbf{G} \\\\ \\midrule \\midrule\n');

% Loop over configurations
for i = 1:height(tableData)
    fprintf(fid, '        \\multirow{4}{*}{%s} ', names(i));
    
    for j = 1:4
        x = round(tableData.Normalised_positions_x(i, j), 3);
        y = round(tableData.Normalised_positions_y(i, j), 3);
        phase = round(tableData.Phases(i, j), 0); % Keep as integer
        mod_eff = round(tableData.Modulation_efficiency(i), 3);
        iwa = round(tableData.IWA(i), 3);
        fwhm = round(tableData.FWHM(i), 3);
        ratio = log10(tableData.Ratio(i));
        
        if j == 1
            fprintf(fid, '& %d & %.3f & %.3f & %d & \\multirow{4}{*}{%.3f} & \\multirow{4}{*}{%.3f} & \\multirow{4}{*}{%.3f} & \\multirow{4}{*}{$10^{%.2f}$} \\\\ \n', j, x, y, phase, mod_eff, iwa, fwhm, ratio);
        else
            fprintf(fid, '        & %d & %.3f & %.3f & %d \\\\ \n', j, x, y, phase);
        end
    end
    if i < height(tableData)
        fprintf(fid, '        \\midrule\n');
    end
end

% Close LaTeX table
fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '    \\caption{Comparison table for five example configurations as presented in Lay \cite{lay_imaging_2005} with normalised results for positions, modulation efficiency (higher is better), inner working angles (lower is better) and resolution angles (lower is better). Ap. is the aperture, $\varphi$ is the assigned phase and G is the nulling ratio (higher exponent is better).}\n');
fprintf(fid, '    \\label{tab:modelling:comparison}\n');
fprintf(fid, '\\end{table}\n');

% Close file
fclose(fid);

end
