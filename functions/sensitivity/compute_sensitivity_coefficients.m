function [C_A, C_phi, C_x, C_y, C_AA, C_Aphi, C_phiphi] = ...
    compute_sensitivity_coefficients(A, phi, B, dB_dx, dB_dy)

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

end