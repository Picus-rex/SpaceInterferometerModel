function export_systematic_coefficients_table(tableData, filename)
% Open file for writing
fid = fopen(filename + ".tex", 'w');

% Write the LaTeX table preamble
fprintf(fid, '    \\begin{tabular}{l*{16}{c}}\n');
fprintf(fid, '        \\toprule\n');
fprintf(fid, '        & \\multicolumn{4}{c}{\\textbf{\\nth{1} order}} & \\multicolumn{12}{c}{\\textbf{\\nth{2} order}} \\\\ \n');
fprintf(fid, '        \\cmidrule(lr){2-5} \\cmidrule(lr){6-17}\n');
fprintf(fid, '        \\textbf{Ap.} & \\textbf{Ap.} & $\\boldsymbol{\\varphi}$ & \\textbf{Pos. x} & \\textbf{Pos. y} & \\multicolumn{4}{c}{\\textbf{Apertures}} & \\multicolumn{4}{c}{\\textbf{Aperture and phase}} & \\multicolumn{4}{c}{\\textbf{Phases}} \\\\ \\midrule \\midrule\n');
fprintf(fid, '        \\multicolumn{4}{l}{\\textit{Star}} & Ap. & 1 & 2 & 3 & 4 & 1 & 2 & 3 & 4 & 1 & 2 & 3 & 4 \\\\ \\cmidrule(lr){1-5} \\cmidrule(lr){6-9} \\cmidrule(lr){10-13} \\cmidrule(lr){14-17}\n');


for i = 1:4
    fprintf(fid, '        %d & %g & %g & %g & %g ', i, round_value(tableData.C_A(i)), round_value(tableData.C_phi(i)), round_value(tableData.C_x(i)), round_value(tableData.C_y(i)));
    
    % First-order coefficients
    for j = 1:4
        fprintf(fid, '& %g ', round_value(tableData.C_AA(i, j)));
    end
    for j = 1:4
        fprintf(fid, '& %g ', round_value(tableData.C_Aphi(i, j)));
    end
    for j = 1:4
        fprintf(fid, '& %g ', round_value(tableData.C_phiphi(i, j)));
    end
    fprintf(fid, '\\\\ \n');
end

fprintf(fid, '        \\midrule\n');
fprintf(fid, '        \\multicolumn{4}{l}{\\textit{Exozodiacal disc}} & Ap. & 1 & 2 & 3 & 4 & 1 & 2 & 3 & 4 & 1 & 2 & 3 & 4 \\\\ \\cmidrule(lr){1-5} \\cmidrule(lr){6-9} \\cmidrule(lr){10-13} \\cmidrule(lr){14-17}\n');

for i = 1:4
    fprintf(fid, '        %d & %g & %g & %g & %g ', i, round_value(tableData.C_A_EZ(i)), round_value(tableData.C_phi_EZ(i)), round_value(tableData.C_x_EZ(i)), round_value(tableData.C_y_EZ(i)));
    
    % Second-order coefficients for Exozodiacal disc
    for j = 1:4
        fprintf(fid, '& %g ', round_value(tableData.C_AA_EZ(i, j)));
    end
    for j = 1:4
        fprintf(fid, '& %g ', round_value(tableData.C_Aphi_EZ(i, j)));
    end
    for j = 1:4
        fprintf(fid, '& %g ', round_value(tableData.C_phiphi_EZ(i, j)));
    end
    fprintf(fid, '\\\\ \n');
end

fprintf(fid, '        \\midrule\n');
fprintf(fid, '        \\multicolumn{4}{l}{\\textit{Local zodiacal disc}} \\\\ \\midrule\n');

for i = 1:4
    fprintf(fid, '        %d & %g ', i, round_value(tableData.C_A_LZ(i)));
    fprintf(fid, '\\\\ \n');
end

fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');

% Close file
fclose(fid);
end





function val = round_value(x)
    if abs(x) > 10
        val = round(x);
    else
        val = round(x, 2);
    end
    if val == 0
        val = 0; % Ensures -0 is displayed as 0
    end
end
