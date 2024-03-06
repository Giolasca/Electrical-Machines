%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ----------------    Transformer Single Phase   ------------------ %%%%
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
% clc; clear; close all;

% Adds the 'functions' directory to the search path
addpath('functions_single_phase');


%% Transformer Specifications
Vp = 10*10^3;       % Primary voltage [V]
Vs = 400;           % Secondary voltage [V]
f = 50;             % Frequency [Hz]
Sn = 1100 * 10^3;   % Nominal apparent power [VA]
mu_0 = 4*pi*1e-7;   % Vacuum permeability


%% Constraints and Design choices
k_thickness = 0.01;     % Thickness of insulation layer [m]
d_core_wire = 0.01;     % Space between the windings and core [m]
air_gap = 0.1e-3;       % Air gap of insulation layer [m]

% Limb type transformer 

% Hypothetical height for the transformer windings [m]
h_windings = 1;  

% Hypothetical thikness for the transformer box [m]
d_tank = 0.001;  


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
density_copper = 8960;       % Density kg/m³
rho_cu = 1.68e-8;      % Resistivity of copper wire (ohm*m)
rho_al = 2.82e-8;      % Resistivity of aluminum wire (ohm*m)

% Choose the material
disp('Choose the material for the core: (1) M400_50  (2) M1000_100');
mat = input('Your choice: ');

switch mat
    case 1
        material = 'M400_50';
        density_steel = density_M400_50;
        cost_steel = 4;
    case 2
        material = 'M1000_100';
        density_steel = density_M1000_100;
        cost_steel = 12;
    otherwise
        disp('Wrong selection')
end

% Choose the efficiency and the rated output
disp('Choose efficiency with rated output power at unity power factor and rated output voltage.');
disp('(1): Efficiency > 85% with the 50% rated output');
disp('(2): Efficiency > 95% with the 75% rated output');
disp('(3): Efficiency > 92% with the 100% rated output');
eff = input('Your choice: ');

switch eff
    case 1
        eff = 0.85;
        load_factor = 0.5;
    case 2
        eff = 0.95;
        load_factor = 0.5;
    case 3
        eff = 0.92;
        load_factor = 0.5;
    otherwise
        disp('Wrong selection')
end


%% Transformer development

disp('Enter the values ​​of B_core, J_winding and K_insulation');
B_core = input('B_core: ');              % Magnetic flux density [T]
J_winding = input('J_windings: ');       % Density current [A*m^2]
K_insulation = input('K_insulation: ');  % Tension insulation [V]

% Trasformation rate
ks = Vp/Vs;
Sn_rated = Sn * load_factor;  % Rated power [VA]

% Calculate primary and secondary current [A]
Ip = Sn_rated/Vp;
Is = Sn_rated/Vs;

% Number of primary and secondary turns
N1 = ceil(Vp/K_insulation);
N2 = ceil(N1 * Vs/Vp);

% Cross-sectional area of wire for primary coil [m^2]
A1_wire = Ip / J_winding;

% Diameter of the core from the cross-sectional area [m]
d1_wire = sqrt(4 * A1_wire / pi);

% Cross-sectional area of wire for secondary coil [m^2]
A2_wire = Is / J_winding;

% Diameter of the core from the cross-sectional area [m]
d2_wire = sqrt(4 * A2_wire / pi);

% Cross-sectional area of the core based on flux density [m^2]
A_core = (Vp / (sqrt(2) * pi * N1 * f * B_core)) * 100/90;

% Diameter of the core from the cross-sectional area [m]
d_core = sqrt(4 * A_core / pi);

% -> Assuming circular winding 
% Length of wire for primary and secondary winding [m]
l1_wire = ceil(N1 * pi * d_core);
l2_wire = ceil(N2 * pi * d_core);



%% Core and winding geometry

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
w_core = 2 *w1 + 2*d_core + 2*w2 + 8*d_core_wire;

% Volume of the core [m^3], 90% fill factor
Vol_core = 2 * A_core * (h_core + (w_core - 2*d_core));

% Core weight [kg]
Kg_core = Vol_core * density_M400_50;

% Radius of the core [m]
R_core = d_core / 2;

% Radii of primary windings [m]
r11_wire = R_core + k_thickness;      % Inner radius primary winding
r12_wire = r11_wire + w1;     % Outer radius primary winding
r1m_wire = r11_wire + w1/2;   % Midpoint radius primary winding

% Radii of secondary windings [m]
rmg = r12_wire + k_thickness/2;       % Midpoint radius between windings
r21_wire = r12_wire + k_thickness;    % Inner radius secondary winding
r22_wire = r21_wire + w2;     % Outer radius secondary winding
r2m_wire = r22_wire + w2/2;   % Midpoint radius secondary winding



%% Resistance of the Core

if (mat == 1)
    Bmax = Bmax_M400_50;
    Ploss = Ploss_M400_50;
    B = B_M400_50;
    H = H_M400_50;
else
    Bmax = Bmax_M1000_100;
    Ploss = Ploss_M1000_100;
    B = B_M1000_100;
    H = H_M1000_100;
end

% Find index in array B_max
[~, index_p] = min(abs(Bmax - B_core)); 

% Steel core power loss [W]
P_loss_steel = Kg_core * Ploss(index_p); 

% Additional power loss in steel (empirically) [W]
P_add = 0.2 * P_loss_steel;

% Find index in array B
[~, index] = min(abs(B - B_core)); 

% Permeability of steel from B/H [Wb/A]
mu_steel = B(index) / H(index); 

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
L12_windings = 2*pi*N1^2*mu_0/h_windings * (r1m_wire*w1/3 + r2m_wire*w2/3 + rmg * k_thickness + w1^2/12 - w2^2/12);

% Reactances for primary and secondary winding [Ohm]
X1 = 2 * pi * f * L1_windings;  
X2 = 2 * pi * f * L2_windings;  
X12 = 2 * pi * f * L12_windings; 

% Impedences of primary and secondary windings [Ohm]
Z1 = R1 + 1i*(X1 + X12);
Z2 = R2 + 1i*(X2 + X12);



%% Efficiency test

% Output power [W]
Pout = Sn_rated;   
            
% Copper losses [W]
Pcu = (Ip^2 * R1 + Is^2 * R2); 

% Input power [W]
Pin = Pout + P_loss_steel + Pcu;  

% Efficiency
eta = Pout / Pin; 

% Display efficiency
disp('Efficiency:');
disp(eta);



%% Costs

% Copper wire cost [€/kg]
Cu_cost = 20;

% Cost contribution oil cooled transformer [€/kg]
Oil_cost = 5;

% Cost Steel [€]
Cost_steel = Kg_core * cost_steel;

% Total volumes of windings [m^3]
Vol_w1 = A1_wire * l1_wire * density_copper;
Vol_w2 = A2_wire * l2_wire * density_copper;

% Cost Copper [€]
Cost_copper1 = Vol_w1 * Cu_cost;
Cost_copper2 = Vol_w2 * Cu_cost;

% Total transformer height [m]
tot_height = h_core;

% Total transformer lenght [m]
tot_lenght = w_core + w1 + w2 + 4*k_thickness;

% Total transformer width [m]
if w1 > w2
    tot_width = 2*w1 + d_core + 2*k_thickness;
else
    tot_width = 2*w2 + d_core + 2*k_thickness;
end

% Total volume of transformer [m^3]
Vol_tot = tot_width * tot_lenght * tot_height; 

% Total volume of oil tank [m^3]
Vol_oil = Vol_tot - Vol_core - Vol_w1 - Vol_w2;

% Cost Oil [€]
Cost_oil = Vol_oil * Oil_cost;

% Cost tank [€]
Cost_tank = 1500;

% Total Cost [€]
Cost_tot = Cost_steel + Cost_copper1 + Cost_copper2 + Cost_oil + Cost_tank;

% Display cost
disp('Cost:');
disp(Cost_tot);