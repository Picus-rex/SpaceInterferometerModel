function [ratio, rejection] = ...
    compute_nulling_ratio(N, phases, positions, theta_star, lambda)

% Defrère D. Characterizing extrasolar planetary systems 
% using infrared interferometry and Absil O. Astrophysical studies of 
% extrasolar planetary systems using infrared interferometric techniques 
% [Internet]. 2006. Available from: 
% https://theses.hal.science/tel-00124720v1/file/Absil06_thesis.pdf



sum = N;  
for j = 1:N
    for k = 1:N        
        if j ~= k
            Delta_x = positions(j,1) - positions(k,1);
            Delta_y = positions(j,2) - positions(k,2);
            dcphi = cos(phases(j) - phases(k));
            b_jk = sqrt(Delta_x^2 + Delta_y^2);
            arg = 2*pi * b_jk * theta_star / lambda;
    
            sum = sum + dcphi * (2 * besselj(1, arg) / arg);
        end
    end
end

ratio = sum / N;

rejection = 1/ratio;




end