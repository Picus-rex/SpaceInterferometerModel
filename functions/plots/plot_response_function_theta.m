function h = plot_response_function_theta(theta_range, RFs)
%PLOT_RESPONSE_FUNCTION_THETA Plots the response functions over a range of
%angular coordinates.
%
% INPUTS:
%   theta_range[Np x 1] Vector of angular coordinates [rad]
%   RFs[Ns x Np]        Matrix of response functions, where each row 
%                       corresponds to a different simulation
%                      
% OUTPUTS:
%   h[handle]           Handle to the plot
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure;
hold on;

for i = 1:size(RFs, 1)
    plot(rad2mas(theta_range), RFs(i, :), "LineWidth", 1.5);
end

xlabel("Observation angle [mas]")
ylabel("Response function")
grid minor;

end