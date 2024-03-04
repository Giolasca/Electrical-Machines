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

% --- CONSTRAINTS↑ AND MATERIAL PROPERTIES↓ END STEP 3 ---
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


% --- END STEP 1 ---
%% Transformer development

% Trasformation rate
ks = Vp/Vs;

% Calculate primary and secondary current [A]
Ip = Sn/Vp;
Is = Sn/Vs;

% Number of primary and secondary turns
N1 = ceil(Vp/Max_V_wire);
N2 = ceil(N1 * Vs/Vp);

% --- END STEP 2 ---

B_core = 1.2;         % Initial guess for magnetic flux density [T]
J_winding = 2.5e6;    % Initial guess for density current [A*m^2]
K_insulation = 0.25;  % Initial guess for tension insulation [V]

% Cross-sectional area of wire for secondary coil [m^2]
A1_wire = Is / J_winding;

% Diameter of the core from the cross-sectional area [m]
d1_wire = sqrt(4 * A1_wire / pi);

% Cross-sectional area of wire for secondary coil [m^2]
A2_wire = Ip / J_winding;

% Diameter of the core from the cross-sectional area [m]
d2_wire = sqrt(4 * A2_wire / pi);

% Cross-sectional area of the core based on flux density [m^2]
A_core = Vp / (sqrt(2) * pi * N1 * f * B_core);

% Diameter of the core from the cross-sectional area [m]
d_core = sqrt(4 * A_core / pi);

% Peak flux density in the core
phi_max = B_core * A_core;

% Assuming circular winding 
% Length of wire for primary and secondary winding [m]
l1_wire = ceil(N1 * pi * d_core);
l2_wire = ceil(N2 * pi * d_core);


% --- END STEP 3&5? ---
%% Core and winding geometry

% Hypothetical height for the transformer windings [m]
h_windings = 1.5;  

% Height of the transformer [m]
h_core = h_windings + 2*d_core;  

% Number of overlpas on the primary
N1_overlaps = round((N1 * (d1_wire + k_thickness)) / h_windings);

% Total final width of the primary winding [m]
w1_winding = N1_overlaps * (d1_wire + k_thickness);

% Number of overlpas on the secondary
N2_overlaps = round((N2 * (d2_wire + k_thickness)) / h_windings);

% Total final width of the secondary winding [m]
w2_winding = N2_overlaps * (d2_wire + k_thickness);

% Total final width of the core [m]
w_core = d_core*2 + 10*k_thickness + 2*w1_winding + 2*w2_winding;

% Volume of the core [m^3], 90% fill factor
V_core = 2 * A_core * (h_core + (w_core - 2*d_core)) * 0.90;

% Core weight [kg]
Kg_core = V_core * density_M400_50;

% Radius of the core [m]
R_core = d_core / 2;

% Radii of primary windings [m]
r11_wire = R_core + k_thickness;      % Inner radius primary winding
r12_wire = r11_wire + w1_winding;     % Outer radius primary winding
r1m_wire = r11_wire + w1_winding/2;   % Midpoint radius primary winding

% Radii of secondary windings [m]
rmg = r12_wire + k_thickness/2;       % Midpoint radius between windings
r21_wire = r12_wire + k_thickness;    % Inner radius secondary winding
r22_wire = r21_wire + w2_winding;     % Outer radius secondary winding
r2m_wire = r22_wire + w2_winding/2;   % Midpoint radius secondary winding


% --- END STEP 4 ---
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

% Permeability of the core material [H/m]
mu = mu_0 * mu_steel;   

% Magnetizing Reactance X0 [Ohm]
X0 = 2 * pi * f * N1^2 * (mu * A_core / l_m); 

% Core impedance [Ohm]
Z0 = (1i * R0 * X0) / (R0 + 1i * X0);



%% Resistances and Reactances of the windings

% Resistance for primary and secondary winding [Ohm]
R1 = rho_cu * l1_wire / A1_wire;  
R2 = rho_cu * l2_wire / A2_wire; 

% Inductances of the windings [H]
L1_windings = N1^2 * A1_wire * mu_steel;  
L2_windings = N2^2 * A2_wire * mu_steel;  
L12_windings = 2*pi*N1^2*mu_0/h_windings * (r1m_wire*w1_winding/3 + r2m_wire*w2_winding/3 + rmg * k_thickness + w1_winding^2/12 - w2_winding^2/12);

% Reactances for primary and secondary winding [Ohm]
X1 = 2 * pi * f * L1_windings;  
X2 = 2 * pi * f * L2_windings;  
X12 = 2 * pi * f * L12_windings; 

% Impedences of primary and secondary windings [Ohm]
Z1 = R1 + 1i*(X1 + X12);
Z2 = R2 + 1i*(X2 + X12);

% disp other variables (dimensions)


% --- END STEP 6 ---
%% Efficiency test

% Load factors
load_factors = [0.5, 0.75, 1]; 
eta = zeros(size(load_factors)); 

for i = 1:length(load_factors)
    Pout = Sn * load_factors(i);     % Output power [W]
    Pcu = (Ip^2 * R1 + Is^2 * R2) * load_factors(i)^2;  % Copper losses [W]
    Pin = Pout + P_loss_steel + Pcu;    % Input power [W]
    eta(i) = Pout / Pin;    % Calculate efficiency
end

% Display efficiencies
disp('Efficiencies for various loads:');
disp(eta);


% --- END STEP 7 ---


%% Costs

% Copper and Aluminum wire cost [€/kg]
Cu_cost = 20;
Al_cost = 4;

% M1000 and M400 steel cost [€/kg]
M_1000 = 4;
M_400 = 12;

% Cost contribution oil cooled transformer [€/kg]
Oil_cost = 5;
tank_cost = 1;

% Calculate the cost of copper or aluminum wire for primary and secondary windings
cost_wire_primary = (rho_cu * l1_wire * Cu_cost);
cost_wire_secondary = (rho_cu * l2_wire * Cu_cost);

% Calculate the cost of steel core
cost_core = Kg_core * M_400;

% Calculate the total cost of the transformer
total_cost = cost_wire_primary + cost_wire_secondary + cost_core;

disp('Total cost of the transformer in €:');
disp(total_cost);


%% Sensitivity Analysis

% Define initial design parameters
Bcore_new = B_core;
Jwinding_new = J_winding;
Kinsulation_new = K_insulation;

% Define sensitivity increments
delta_Bcore = 0.1;  % Change in peak flux density [T]
delta_Jwinding = 0.1e6;  % Change in winding current density [A/m^2]
delta_Kinsulation = 0.1;  % Change in insulation voltage per turn [V]

% Maximum number of iterations
max_iterations = 100;
iteration = 0;

% Initialize flag for optimal solution
optimal_solution_found = false;

while ~optimal_solution_found && iteration < max_iterations
    % Recalculate the total cost with new design parameters
    % (repeat the calculations for wire cost, core cost, and total cost)
    cost_wire_primary_new = (rho_cu * l1_wire * Cu_cost) / 1000; % Convert kg to g
    cost_wire_secondary_new = (rho_cu * l2_wire * Cu_cost) / 1000; % Convert kg to g
    cost_core_new = Kg_core * M_400;
    total_cost_new = cost_wire_primary_new + cost_wire_secondary_new + cost_core_new;

    % Compare the total cost with the original cost and analyze the impact
    if total_cost_new < total_cost
        optimal_solution_found = true;
    else
        % Update design parameters for the next iteration
        Bcore_new = Bcore_new + delta_Bcore;
        Jwinding_new = Jwinding_new + delta_Jwinding;
        Kinsulation_new = Kinsulation_new + delta_Kinsulation;
    end
    
    % Increment iteration counter
    iteration = iteration + 1;
end

% Display message if optimal solution is found
if optimal_solution_found
    disp('New design parameters have lowered the total cost.');
    disp('Adopting these parameters for optimal cost.');
    % Display the final parameters
    disp('FINAL PARAMETERS:');
    disp(['Peak Flux Density (Bcore_new): ', num2str(Bcore_new), ' T']);
    disp(['Winding Current Density (Jwinding_new): ', num2str(Jwinding_new), ' A/m^2']);
    disp(['Insulation Voltage per Turn (Kinsulation_new): ', num2str(Kinsulation_new), ' V']);
else
    disp('Maximum number of iterations reached. No new optimal solution found.');
end

%% Three-Phase Integrated Design (Y-connected Transformer)

% Total inductance
L_leakage = 2 * L12_windings; % Total leakage inductance

% Update core parameters for three-phase design
w_core_3phase = w_core;  % Total width of the three-phase core remains the same
h_core_3phase = h_core;  % Height remains the same
d_core_3phase = d_core;  % Depth remains the same

% Update winding parameters for three-phase design
N1_3phase = N1 * sqrt(3);  % Total number of primary turns in three-phase (Y connection)
N2_3phase = N2 * sqrt(3);  % Total number of secondary turns in three-phase (Y connection)

% Calculate leakage impedance per phase
% Assuming equal leakage impedance for each phase
Z_leakage_per_phase = 1i * 2 * pi * f * L_leakage; % Assuming Y-connected transformer

% Update leakage impedance for three-phase
Z_leakage_3phase = Z_leakage_per_phase / 3;

% Update magnetizing branch parameters for three-phase
R0_3phase = R0;  % Magnetization resistance for three-phase remains the same
X0_3phase = X0;  % Magnetization reactance for three-phase remains the same

% Update total cost for three-phase design
cost_core_3phase = cost_core;                      % Total core cost for three-phase remains the same
cost_wire_primary_3phase = cost_wire_primary;      % Total primary wire cost for three-phase remains the same
cost_wire_secondary_3phase = cost_wire_secondary;  % Total secondary wire cost for three-phase remains the same
total_cost_3phase = total_cost;                    % Total cost for three-phase design remains the same

% Display the updated parameters for three-phase design
disp('Updated Parameters for Three-Phase Integrated Design:');
disp(['Total Width of Core: ', num2str(w_core_3phase)]);
disp(['Total Height of Core: ', num2str(h_core_3phase)]);
disp(['Total Depth of Core: ', num2str(d_core_3phase)]);
disp(['Total Number of Primary Turns: ', num2str(N1_3phase)]);
disp(['Total Number of Secondary Turns: ', num2str(N2_3phase)]);
disp(['Leakage Impedance per Phase: ', num2str(Z_leakage_per_phase)]);
disp(['Magnetization Resistance for Three-Phase: ', num2str(R0_3phase)]);
disp(['Magnetization Reactance for Three-Phase: ', num2str(X0_3phase)]);
disp(['Total Cost for Three-Phase Design: ', num2str(total_cost_3phase)]);
