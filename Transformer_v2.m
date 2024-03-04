%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ----------------------    Transformer    ------------------------ %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimized design of a three phase YD distribution transformer 
% with nominal primary to secondary line to line voltages 10kV/400V,
% 50Hz with nominal apparent power 1100kVA. 

% Efficiency > 85% with the 50% rated output power 
% at unity power factor and rated output voltage.

% Efficiency > 95% with the 75% rated output power 
% at unity power factor and rated output voltage.

% Efficiency > 92% with the 100% rated output power 
% at unity power factor and rated output voltage.


%% Initialization
clc; clear; close all;
addpath('functions');


%% Constants
mu_0 = 4*pi*1e-7;   % Vacuum permeability


%% Transformer Specifications
Vp = 10*10^3;       % Primary voltage [V]
Vs = 400;           % Secondary voltage [V]
f = 50;             % Frequency [Hz]
Sn = 1100 * 10^3;   % Nominal apparent power [VA]


%% Constraints
J_Cu_air = 1;    % Current density copper [A/mm^2] (rms)
J_Cu_oil = 2;    % Current density copper [A/mm^2] (rms)

Max_V_wire = 3;          % Max voltage handle by wire insulation [V]
k_thickness = 0.01;     % Thickness of insulation layer [m]

air_gap = 0.1e-3;       % Air gap of insulation layer [m]

%% Material properties

% M400-50 Alloy 
file_M400_50 = 'M400-50_data.xlsx';   % Set the path of the Excel file
M400_50 = readtable(file_M400_50);

% Extract the columns as separate cells
H_M400_50 = M400_50{:, 1};           % Extract H
B_M400_50 = M400_50{:, 2};           % Extract B
f_M400_50 = M400_50{1:18, 4};        % Extract f
Bmax_M400_50 = M400_50{1:18, 5};     % Extract Bmax
Ploss_M400_50 = M400_50{1:18, 6};    % Extract Ploss

density_M400_50 = 77000;          % Density kg/m³


% M1000-100 Alloy 
file_M1000_100 = 'M1000-100_data.xlsx';   % Set the path of the Excel file
M1000_100 = readtable(file_M1000_100);

% Extract the columns as separate cells
H_M1000_100 = M1000_100{:, 1};           % Extract H
B_M1000_100 = M1000_100{:, 2};           % Extract B
f_M1000_100 = M1000_100{1:18, 4};        % Extract f
Bmax_M1000_100 = M1000_100{1:18, 5};     % Extract Bmax
Ploss_M1000_100 = M1000_100{1:18, 6};    % Extract Ploss

density_M1000_100 = 78000;        % Density kg/m³


% The resistivity of copper or aluminum depends on the temperature 
% conditions and purity of the copper used in the transformer

% Wire Material
rho_cu = 1.68e-8;      % Resistivity of copper wire (ohm*m)
rho_al = 2.82e-8;      % Resistivity of aluminum wire (ohm*m)



%% Transformer development

% Trasformation rate
ks = Vp/Vs;

% Calculate primary and secondary current [A]
Ip = Sn/Vp;
Is = Sn/Vs;

% Number of primary and secondary turns
N1 = ceil(Vp/Max_V_wire);
N2 = ceil(N1 * Vs/Vp);

B_core = 1.2;         % Initial guess for magnetic flux density [T]
J_winding = 2.5e6;    % Initial guess for density current [A*m^2]
K_insulation = 0.25;  % Initial guess for tension insulation [V]



%% Core and winding geometry

% Cross-sectional area [m^2] and diameter [m] of core based on flux density 
[A_core, d_core] = geom_core(Vp, N1, f, B_core);

% Cross-sectional area [m^2], diameter [m] and length [m] of wires 
[A1_wire, d1_wire, l1_wire] = geom_wire(Is, J_winding, N1, d_core);
[A2_wire, d2_wire, l2_wire] = geom_wire(Ip, J_winding, N2, d_core);


% Hypothetical height for the transformer windings [m]
h_windings = 1.5;  

% Height of the transformer [m]
h_core = h_windings + 2*d_core;  

% Number of overlpas on the primary
N1_overlaps = round((N1 * (d1_wire + k_thickness)) / h_windings);

% Total final width of the primary winding [m]
w1 = N1_overlaps * (d1_wire + k_thickness);

% Number of overlpas on the secondary
N2_overlaps = round((N2 * (d2_wire + k_thickness)) / h_windings);

% Total final width of the secondary winding [m]
w2 = N2_overlaps * (d2_wire + k_thickness);

% Total final width of the core [m]
w_core = d_core*2 + 10*k_thickness + 2*w1 + 2*w2;

% Volume of the core [m^3], 90% fill factor
V_core = 2 * A_core * (h_core + (w_core - 2*d_core)) * 0.90;

% Core weight [kg]
Kg_core = V_core * density_M400_50;





%% Resistance of the Core

% Find index in array B_max
[~, index_p] = min(abs(Bmax_M400_50 - B_core)); 

% Steel core power loss [W]
P_loss_steel = Kg_core * Ploss_M400_50(index_p); 

% Additional power loss in steel (empirically) [W]
P_add = 0.2 * P_loss_steel;

% Find index in array B
[~, index] = min(abs(B_M400_50 - B_core)); 

% Permeability of steel from B/H [Wb/A]
mu_steel = B_M400_50(index) / H_M400_50(index); 

% Magnetic circuit length [m]
l_m = (w_core - d_core) * 2 + (h_core - d_core) * 2; 

% Magnetizing current [A]
I_M = (B_core * l_m) / ((N1 + N2) * mu_steel); 

% Magnetization resistance R0 [Ohm]
R0 = Vp^2 / P_loss_steel; 



% %% Open Circuit Test (No Load)
% 
% % Power loss in the copper [W]
% P_cu_M = R1 * I_M^2;
% 
% % Total active power loss [W]
% P0 = P_cu_M + P_add + P_loss_steel;
% 
% % Total reactive power [VAr]
% Q0 = sqrt(Sn^2 - P0^2);
% 
% % Core reactance [Ohm]
% X0 = Vp^2 / Q0;
% 
% % Core impedance [Ohm]
% Z0 = (1i * R0 * X0) / (R0 + 1i * X0);



%% Resistances and Reactances of the windings

% Resistance [Ohm]  -  Inductance [H]  -  Reactance [Ohm]
[R1, L1, X1] = winding(rho_cu, A1_wire, l1_wire, N1, mu_steel, f);
[R2, L2, X2] = winding(rho_cu, A1_wire, l1_wire, N1, mu_steel, f);
[L12, X12] = mutual_inductance(d_core, k_thickness, h_windings, w1, w2, N1, f);

% Impedences of primary and secondary windings [Ohm]
Z1 = R1 + 1i*(X1 + X12);
Z2 = R2 + 1i*(X2 + X12);



%% Efficiecy test

for i = 1:length(load_factors)
    Pout = Sn * load_factors(i); % Potenza in uscita
    Pcu = (Ip^2 * R1 + Is^2 * R2) * load_factors(i)^2; % Perdite nel rame
    Pin = Pout + P_loss_steel + Pcu; % Potenza in ingresso
    efficiencies(i) = Pout / Pin; % Calcolo efficienza
end

% Stampa delle efficienze
disp('Efficienze per i vari carichi:');
disp(efficiencies);




%% Costs