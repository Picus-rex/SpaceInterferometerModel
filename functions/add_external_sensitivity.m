function epsilon = add_external_sensitivity(instrument, environment)
%ADD_EXTERNAL_SENSITIVITY Consider the effects of environmental leakage on
%the behaviour of a nulling interferometer. This function requires inputs
%structured following the configuration files.
%
% INPUTS:
%   instruments[struct] Fields linked to the interferometer. Includes:
%                        - positions    Nx2
%                        - lambda       1
%                        - diameter     1
%                        - phase_shifts Nx1
%                        - apertures    1
%                        - intensities  Nx1
%   environment[struct] Fields linked to the environmental conditions. For
%                       details, consider the configuration file
%
% OUTPUTS:
%   epsilon[1]          Perturbation associated to the response function.
%
% REFERENCES:
%   Lay OP. Systematic errors in nulling interferometers. 2004;  
%
% NOTES:
%   - This function is part of Monte Carlo analysis.
%
% VERSION HISTORY:
%   2025-02-27 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract data from instrument
positions = instrument.positions;
lambda = instrument.wavelength;
diameter = instrument.diameter(1);
phi = instrument.phase_shifts;
N = instrument.apertures;
A = instrument.intensities';

% Extract data from perturbations
Fstar = environment.disturbances.star_flux / diameter;       
theta_star = environment.stellar_angular_radius;             

Fplanet = environment.disturbances.planet_flux; 
theta_planet_x = cell2mat(environment.exoplanet_position(1));
theta_planet_y = cell2mat(environment.exoplanet_position(2));

FEZ = environment.disturbances.exozodiacal_flux; 
theta_ez = environment.disturbances.exozodiacal_extension;

BLZ = environment.disturbances.localzodiacal_flux;
DeltaOmega = environment.disturbances.effective_solid_angle;

% Perturbations
perturbations = environment.disturbances.perturbations;
delta_a   = perturbations.intensity * randn(N,1);      
delta_phi = perturbations.phase * randn(N,1);      
delta_x   = perturbations.x_position * randn(N,1);     
delta_y   = perturbations.y_position * randn(N,1);    

% Star Contributions
B_star = zeros(N, N);
dB_dx = zeros(N, N);
dB_dy = zeros(N, N);

for j = 1:N
    for k = 1:N
        xjk = positions(j,1) - positions(k,1);
        yjk = positions(j,2) - positions(k,2);
        b_jk = sqrt(xjk^2 + yjk^2);
        arg = 2*pi * b_jk * theta_star / lambda;
        
        if j ~= k
            B_star(j,k) = 2 * Fstar * besselj(1, arg) / arg;
            dB_dx(j,k) = -2 * xjk * Fstar * besselj(2, arg) / b_jk^2;
            dB_dy(j,k) = -2 * yjk * Fstar * besselj(2, arg) / b_jk^2;
        else
            B_star(j,k) = 2 * Fstar * 0.5;
            dB_dx(j,k) = -2 * xjk * Fstar * 0;
            dB_dy(j,k) = -2 * yjk * Fstar * 0;
        end
        
    end
end

[C_A, C_phi, C_x, C_y, C_AA, C_Aphi, C_phiphi] = ...
    compute_coefficients(A, phi, B_star, dB_dx, dB_dy, diameter);

N_star_ideal = 0;
for j = 1:N
    for k = 1:N
        N_star_ideal = N_star_ideal + ...
            A(j) * A(k) * cos(phi(j)-phi(k)) * B_star(j,k);
    end
end

% Resulting perturbations
deltaN_star_1 = sum(C_A .* delta_a) + sum(C_phi .* delta_phi) + ...
    sum(C_x .* delta_x) + sum(C_y .* delta_y);
deltaN_star_2 = 0;
for j = 1:N
    for k = 1:N
        deltaN_star_2 = deltaN_star_2 + C_AA(j,k)*delta_a(j)*delta_a(k) ...
            + C_Aphi(j,k)*delta_a(j)*delta_phi(k) + ...
            C_phiphi(j,k)*delta_phi(j)*delta_phi(k);
    end
end

deltaN_star = deltaN_star_1 + deltaN_star_2;

% Exozodiacal Contributions

B_EZ = zeros(N, N);
dB_dx = zeros(N, N);
dB_dy = zeros(N, N);

for j = 1:N
    for k = 1:N
        xjk = positions(j,1) - positions(k,1);
        yjk = positions(j,2) - positions(k,2);
        b_jk = sqrt(xjk^2 + yjk^2);
        arg = 2*pi * b_jk * theta_ez / lambda;
        
        if j ~= k
            B_EZ(j,k) = 2 * FEZ * besselj(1, arg) / arg;
            dB_dx(j,k) = -2 * xjk * FEZ * besselj(2, arg) / b_jk^2;
            dB_dy(j,k) = -2 * yjk * FEZ * besselj(2, arg) / b_jk^2;
        else
            B_EZ(j,k) = 2 * FEZ * 0.5;
            dB_dx(j,k) = -2 * xjk * FEZ * 0;
            dB_dy(j,k) = -2 * yjk * FEZ * 0;
        end
 
    end
end

[C_A_EZ, C_phi_EZ, C_x_EZ, C_y_EZ, C_AA_EZ, C_Aphi_EZ, C_phiphi_EZ] = ...
    compute_coefficients(A, phi, B_EZ, dB_dx, dB_dy, diameter);

N_EZ_ideal = 0;
for j = 1:N
    for k = 1:N
        N_EZ_ideal = N_EZ_ideal + A(j)*A(k)*cos(phi(j)-phi(k))*B_EZ(j,k);
    end
end 

deltaN_EZ_1 = sum(C_A_EZ .* delta_a) + sum(C_phi_EZ .* delta_phi) + ...
    sum(C_x_EZ .* delta_x) + sum(C_y_EZ .* delta_y);
deltaN_EZ_2 = 0;
for j = 1:N
    for k = 1:N
        deltaN_EZ_2 = deltaN_EZ_2 + C_AA_EZ(j,k)*delta_a(j)*delta_a(k) ...
            + C_Aphi_EZ(j,k)*delta_a(j)*delta_phi(k) + ...
            C_phiphi_EZ(j,k)*delta_phi(j)*delta_phi(k);
    end
end

deltaN_EZ = deltaN_EZ_1 + deltaN_EZ_2;

% Local Zodiacal Contributions

N_LZ_ideal = BLZ * sum(A.^2 .* DeltaOmega);
C_A_LZ = 2 * BLZ * (A.^2 .* DeltaOmega);  
deltaN_LZ = sum(C_A_LZ .* delta_a);

% Planet Contributions

N_planet_ideal = 0;
for j = 1:N
    for k = 1:N
        phase_term = (2*pi/lambda) * ...
            ((positions(j,1)-positions(k,1)) * theta_planet_x + ...
            (positions(j,2)-positions(k,2)) * theta_planet_y);

        N_planet_ideal = N_planet_ideal + A(j) * A(k) * ...
            (cos(phi(j)-phi(k)) * cos(phase_term) - ...
            sin(phi(j)-phi(k)) * sin(phase_term));
    end
end
N_planet_ideal = Fplanet * N_planet_ideal;

% Total

N0 = N_star_ideal + N_EZ_ideal + N_LZ_ideal + N_planet_ideal;
deltaN_total = deltaN_star + deltaN_EZ + deltaN_LZ;

% Fractional (relative) perturbation
epsilon = deltaN_total / N0;


end





function [C_A, C_phi, C_x, C_y, C_AA, C_Aphi, C_phiphi] = ...
    compute_coefficients(A, phi, B, dB_dx, dB_dy, diameter)

N = length(A);

% 1st order
C_A = zeros(N,1);
C_phi = zeros(N,1);
C_x = zeros(N,1);
C_y = zeros(N,1);

for j = 1:N

    sumA = 0; 
    sumPhi = 0; 
    sumX = 0; 
    sumY = 0;
    
    for k = 1:N
        sumA = sumA + A(k) * cos(phi(j)-phi(k)) * B(j,k);
        if j ~= k 
            sumPhi = sumPhi + A(k) * sin(phi(j)-phi(k)) * B(j,k);
        end
        sumX = sumX + A(j)*A(k)*cos(phi(j)-phi(k)) * dB_dx(j,k);
        sumY = sumY + A(j)*A(k)*cos(phi(j)-phi(k)) * dB_dy(j,k);
    end

    C_A(j) = 2 * A(j) * sumA;
    C_phi(j) = -2 * A(j) * sumPhi;
    C_x(j) = 2 * sumX;
    C_y(j) = 2 * sumY;
end

% 2nd order
C_AA = zeros(N,N);
C_Aphi = zeros(N,N);
C_phiphi = zeros(N,N);

for j = 1:N
    for k = 1:N
        
        C_AA(j,k) = A(j) * A(k) * cos(phi(j)-phi(k)) * B(j,k);

        if j ~= k
            
            C_Aphi(j,k) = 2 * A(j) * A(k) * sin(phi(j)-phi(k)) * B(j,k);
            C_phiphi(j,k) = A(j) * A (k) * cos(phi(j)-phi(k)) * B(j,k);

        else

            tempAphi = 0;
            tempPhiPhi = 0;

            for l = 1:N
                if l ~= j
                    tempAphi = tempAphi + A(l)*sin(phi(j)-phi(l)) * B(k,l);
                    tempPhiPhi = tempPhiPhi + A(l)*cos(phi(j)-phi(l))*B(j,l);
                end
            end

            C_Aphi(j,j) = -2 * A(j) * tempAphi;
            C_phiphi(j,j) = -A(j) * tempPhiPhi;
        end
    end
end

C_AA = C_AA * diameter;
C_Aphi = C_Aphi * diameter;
C_phiphi = C_phiphi * diameter;

end