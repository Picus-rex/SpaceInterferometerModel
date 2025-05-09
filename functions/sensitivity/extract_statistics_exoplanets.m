function sumtable = extract_statistics_exoplanets(exotable, name)
%EXTRACT_STATISTICS_EXOPLANETS Evaluate statistic data from exoplanets
%table to study the yields of different planets.
%
% INPUTS:
%   exotable[table]     Exoplanets table from P-Pop with additional fields
%                       from get_ppop_yield. In particular:
%       Nuniverse[Nex1] Universe associated to that planet.
%       Rp[Ne x 1]      Radius of the exoplanet. [R_Earth]
%       Tp[Ne x 1]      Equilibrium temperature of exoplanet. [K]
%       Ds[Ne x 1]      Distance to Sun of host star. [pc]
%       Stype[Ne x 1]   Type of the host star (like 'F', 'G', 'K', 'M').
%       yields[Ne x Ns] Number of detected exoplanet for each simulation
%                       and for each universe of the simulation.
%       name[string]    Name associated to the current array, used as first
%                       column in the output sumtable. 
%
% OUTPUTS:
%   sumtable[table]     Statistic table associated to the given array, with
%                       fields:
%       name[string]    As given by input argument
%       MeanDetections[1]
%       MedianDetections[1]
%       MedianRadius[1]
%       MeanRocky[1]
%       MeanSuperEarth[1]
%       MeanSubNeptune[1]
%       MeanSubJovian[1]
%       MeanJovian[1]
%       MeanPlanets_AFGK[1]
%       MeanPlanets_Mstar[1]
%       MeanTemp[1]
%       MeanDist[1]
%       MinDist[1]
%
% NOTE:
%   Table structure inspired by Dandumont. 
%
% REFERENCES:
%   Kopparapu RK, Hébrard E, Belikov R, Batalha NM, Mulders GD, Stark C, et 
%   al. Exoplanet Classification and Yield Estimates for Direct Imaging 
%   Missions. ApJ. 2018 Apr 1;856(2):122. 
%
% VERSION HISTORY:
%   2025-05-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract elements from tables. Since universes start from 0, add 1 to have
% the exact number of universes.
universes = unique(exotable.Nuniverse);
Nu = max(universes) + 1;  
Ns = size(exotable.yields, 2);  

% Preallocate the results. 
detCounts = zeros(Nu, Ns);
classCounts = struct(...
    'rocky', zeros(Nu, Ns), ...
    'superE', zeros(Nu, Ns), ...
    'subN', zeros(Nu, Ns), ...
    'subJ', zeros(Nu, Ns), ...
    'Jovian', zeros(Nu, Ns), ...
    'AFGK', zeros(Nu, Ns), ...
    'Mstar', zeros(Nu, Ns) ...
);

% For every universe... 
for u = 0:Nu-1

    % Extract the rows associated to the current universe and elements of
    % the planets that are interested.
    idx = exotable.Nuniverse == u;
    Rp_u    = exotable.Rp(idx);
    Tp_u    = exotable.Tp(idx);
    Ds_u    = exotable.Ds(idx);
    Stype_u = exotable.Stype(idx);
    
    % Extraction of seen exoplanets 
    Y = exotable.yields(idx, :);  % Npu x Ns 
    detCounts(u+1, :) = sum(Y, 1);
    
    % Classification by Kopparapu, 2018
    classCounts.rocky(u+1,:)   = sum(Y(Rp_u < 1, :), 1);
    classCounts.superE(u+1,:)  = sum(Y(Rp_u >= 1   & Rp_u < 1.75, :), 1);
    classCounts.subN(u+1,:)    = sum(Y(Rp_u >= 1.75& Rp_u < 3.5,  :), 1);
    classCounts.subJ(u+1,:)    = sum(Y(Rp_u >= 3.5 & Rp_u < 6,    :), 1);
    classCounts.Jovian(u+1,:)  = sum(Y(Rp_u >= 6   & Rp_u < 14.5, :), 1);
    
    % Star classification AFGK vs M
    isAFGK = ismember(Stype_u, {'A','F','G','K'});
    classCounts.AFGK(u+1,:)    = sum(Y(isAFGK, :), 1);
    classCounts.Mstar(u+1,:)   = sum(Y(~isAFGK, :), 1);
    
    % Temperatures and distances for every universes
    for s = 1:Ns
        seen_idx = Y(:,s)==1;
        Tp_seen{u+1, s} = Tp_u(seen_idx);
        Ds_seen{u+1, s} = Ds_u(seen_idx);
    end
end

% Global characteristics, for universe and simulations, by moving to a
% vectorial dimension. 
allDet = detCounts(:);
meanDet = mean(allDet);     
medianDet = median(allDet);

% Average for all classes (except for min)
meanRocky  = mean(classCounts.rocky(:)) ;
meanSuperE = mean(classCounts.superE(:));
meanSubN   = mean(classCounts.subN(:))  ;
meanSubJ   = mean(classCounts.subJ(:))  ;
meanJovian = mean(classCounts.Jovian(:));
meanAFGK   = mean(classCounts.AFGK(:))  ;
meanMstar  = mean(classCounts.Mstar(:)) ;

% Concatenation for Tp and Ds
allTp = vertcat(Tp_seen{:});
allDs = vertcat(Ds_seen{:});
meanTp = mean(allTp);
meanDs = mean(allDs);
minDs  = min(allDs);

% Export
sumtable = table(...
    name, meanDet, medianDet, median(exotable.Rp), ...
    meanRocky, meanSuperE, meanSubN, meanSubJ, meanJovian, ...
    meanAFGK, meanMstar, ...
    meanTp, meanDs, minDs, ...
    'VariableNames', {...
      'Name', 'MeanDetections','MedianDetections','MedianRadius',...
      'MeanRocky','MeanSuperEarth','MeanSubNeptune','MeanSubJovian','MeanJovian',...
      'MeanPlanets_AFGK','MeanPlanets_Mstar',...
      'MeanTemp','MeanDist','MinDist'});

end
