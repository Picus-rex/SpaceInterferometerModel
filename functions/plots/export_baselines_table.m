function export_baselines_table(xArrayTable, linearArrayTable, filename)

% Open file for writing
fid = fopen(filename+".tex", 'w');

% Write the LaTeX table preamble
fprintf(fid, '\\begin{table}[htb!]\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\small\n');
fprintf(fid, '    \\begin{tabular}{clcccc}\n');
fprintf(fid, '        \\toprule\n');
fprintf(fid, '        \\textbf{Distinct baselines} & \\textbf{Type} & \\textbf{Length} [m] & \\textbf{Phase diff.} [deg] & \\textbf{Factor} & \\textbf{Contributions} \\\\\n');
fprintf(fid, '        \\midrule \\midrule\n');

% Write X-Array section
fprintf(fid, '        \\multicolumn{6}{l}{\\textit{X-Array}} \\\\\n');
fprintf(fid, '        \\midrule\n');
write_table_section(fid, xArrayTable);
fprintf(fid, '        \\midrule\n');

% Write Linear Array section
fprintf(fid, '        \\multicolumn{6}{l}{\\textit{Linear Array}} \\\\\n');
fprintf(fid, '        \\midrule\n');
write_table_section(fid, linearArrayTable);

% Close LaTeX table
fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '    \\caption{Classification of baselines for the two proposed arrays in terms of type, length, and contributing repeating baselines for the considered arrays. The numbering associated to the contributing apertures is referenced from the presented configurations in Figure~\ref{fig:modelling:configurations}.}\n');
fprintf(fid, '    \\label{tab:baseline_classifications}\n');
fprintf(fid, '\\end{table}\n');

% Close file
fclose(fid);

end





function write_table_section(fid, tableData)

% Loop over rows in the table and write data to file
for i = 1:height(tableData)
    % Extract relevant data
    length = tableData.B(i);
    phase_diff = rad2deg(tableData.Delta_phi(i));
    factor = tableData.C_i(i);

    idx = find(tableData.contributing(i, :) ~= "");
    contributions = strjoin(tableData.contributing(i, idx), ", ");
    
    if tableData.Imaging_Flag(i) == 0
        type = "Nulling";
    elseif tableData.Imaging_Flag(i) == 1
        type = "Imaging";
    else
        type = "Mixed";
    end


    % Write row to LaTeX file
    fprintf(fid, '        %d & %s & %.3f & %.0f & %.2f & %s \\\\\n', i, type, length, phase_diff, factor, contributions);
end

end
