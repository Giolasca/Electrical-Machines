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
clc; clear; close all;

% Adds the 'functions' directory to the search path
addpath('functions');


%% Constraints and Costants
mu_0 = 4*pi*1e-7;   % Vacuum permeability

k_thickness = 0.0005 * 2;   % Thickness of insulation layer [m]
d_core_wire = 0.01;         % Space between the windings and core [m]
air_gap = 0.1e-3;           % Air gap of insulation layer [m]

% The resistivity of copper or aluminum depends on the temperature 
% conditions and purity of the copper used in the transformer

% Wire Material
density_copper = 8960;    % Density kg/m³
rho_cu = 1.68e-8;         % Resistivity of copper wire (ohm*m)

% Core Material
density_M400_50 = 77000;       % Density kg/m³
density_M1000_100 = 78000;     % Density kg/m³


%% Design choices

% Limb type transformer 

% Hypothetical height for the transformer windings [m]
h_windings = 1;  

% Hypothetical thikness for the transformer box [m]
d_tank = 0.001;  

% Choose the material
disp('Choose the material for the core: (1) M400_50  (2) M1000_100');
mat = input('Your choice: ');

switch mat
    case 1
        material = 'M400_50';
        density_steel = density_M400_50;
    case 2
        material = 'M1000_100';
        density_steel = density_M1000_100;
    otherwise
        disp('Wrong selection')
end


%% Transformer development

Vp = 10*10^3;       % Primary voltage [V]
Vs = 400;           % Secondary voltage [V]
f = 50;             % Frequency [Hz]
Sn = 1100 * 10^3;   % Nominal apparent power [VA]



%% Optimization problem
% Definition of ranges
B_core_range = 0.5:0.1:1.8;             % From 0.5 T to 2 T
J_winding_range = 2.5e6:2.5e6:20e6;   % From 2.5e6 A/m^2 to 20e6 A/m^2
K_insulation_range = 0.25:0.25:3;     % From 0.25 V to 3 V

% Array of required efficiencies and corresponding load factors
Eff_required = [85, 95, 92] / 100;    % Conversion to decimal form
Load_factors = [0.5, 0.75, 1];        % Load factors 50%, 75%, 100%

% Variables to keep track of the best result
min_cost = inf;
best_B_core = 0;
best_J_winding = 0;
best_K_insulation = 0;
best_eta = 0;

% Iteration over every possible combination of B_core, J_winding, and K_insulation
for B_core = B_core_range
    for J_winding = J_winding_range
        for K_insulation = K_insulation_range
            for i = 1:length(Load_factors)
                
                ks = Vp/Vs;    % Trasformation rate
                Sn_rated = Sn * Load_factors(i);  % Rated power [VA]

                % Calculate primary and secondary current [A]
                Ip = Sn/Vp;
                Is = Sn/Vs;

                % Number of primary and secondary coil
                N1 = ceil(Vp/K_insulation);    
                N2 = ceil(N1 * Vs/Vp);         

                % Cross-sectional area [m^2] and diameter [m] of core based on flux density 
                [A_core, d_core] = geom_core(Vp, N1, f, B_core);

                % Cross-sectional area [m^2], diameter [m], and length [m] of wires 
                [A1_wire, d1_wire, l1_wire] = geom_wire(Ip, J_winding, N1, d_core);
                [A2_wire, d2_wire, l2_wire] = geom_wire(Is, J_winding, N2, d_core);

                % Volumes of wire [m^3]
                Vol1_wire = A1_wire * l1_wire;
                Vol2_wire = A2_wire * l2_wire;

                % Weight of the wires [Kg]
                Kg_wire1 = Vol1_wire * density_copper;
                Kg_wire2 = Vol2_wire * density_copper;

                % Final width of the primary and secondary winding [m]
                w1 =  width_windings(N1, d1_wire, k_thickness, h_windings); 
                w2 =  width_windings(N2, d2_wire, k_thickness, h_windings);  
            
                % Height and width of the transformer [m]
                h_core = h_windings + 2*d_core;  
                w_core = 2 *w1 + 2*d_core + 2*w2 + 8*d_core_wire;

                % Volume of the core [m^3] with 90% fill factor and core weight [kg]
                Vol_core = 2 * A_core * (h_core + (w_core - 2*d_core));
                Kg_core = Vol_core * density_steel;
    
                % Magnetic circuit length [m]
                lenght_magnetic = 2 * w_core + 2 * h_core - 4 * d_core; 


                %% Resistances and Reactances of the Core

                % Steel core power loss [W] and Permeability of steel from B/H [Wb/A]
                [Ploss_steel, mu_steel] = steel(material, B_core, Kg_core);

                % Magnetization resistance R0 [Ohm]
                R0 = Vp^2 / Ploss_steel; 

                % Magnetizing Reactance X0 [Ohm]
                X0 = 2 * pi * f * (N1 + N2)^2 * (mu_0 * mu_steel * A_core / lenght_magnetic); 

                % Core impedance (R0 || X0) [Ohm]
                Z0 = (R0 * 1i * X0) / (R0 + 1i * X0);


                %% Resistances and Reactances of the Windings

                % Resistance [Ohm]  -  Inductance [H]  -  Reactance [Ohm]
                [R1, L1, X1] = winding(rho_cu, A1_wire, l1_wire, N1, mu_steel, f);
                [R2, L2, X2] = winding(rho_cu, A2_wire, l2_wire, N2, mu_steel, f);
                [L12, X12] = mutual_inductance(d_core, k_thickness, h_windings, w1, w2, N1, f);

                % Impedences of primary and secondary windings [Ohm]
                Z1 = R1 + 1i*(X1 + X12);
                Z2 = R2 + 1i*(X2 + X12);

                
                %% Short-circuit current referred to the primary side

                % Reverts the secondary impedance to the primary side
                Z2_primary = Z2 * (ks^2);  
                
                % Total short-circuit impedance seen from the primary side
                Z_sc_primary = Z1 + Z2_primary;
                
                % Short-circuit current on the primary side
                I_sc_primary = Vp / abs(Z_sc_primary);

                % Three times the primary nominal current
                Ip_3x = 3 * Ip;

                % Check requirement
                if I_sc_primary > Ip_3x
                    % If it exceeds, this configuration does not meet the requirement
                    continue; % Exclude this configuration
                end


                %% Costs
            
                % Total cost of metals [€]
                [Cost_steel, Cost_copper] = cost_metals(material, Kg_core, Kg_wire1, Kg_wire2);

                % Total volume of transformer [m^3]
                [Vol_tank, Sup_tank] = tank(w_core, h_core, d_core, w1, w2, k_thickness);

                % Total cost of oil and tank [€]
                [Cost_oil, Cost_tank] = cost_oil(Vol_tank, Vol_core, Vol1_wire, Vol2_wire, Sup_tank, d_tank);

                % Total cost of transformer [€]
                if (J_winding > 10e6)
                    Cost_tot = Cost_steel + Cost_copper + Cost_tank + Cost_oil;
                else
                    Cost_tot = Cost_steel + Cost_copper;
                end

            
                %% Efficiency test
            
                % Output power [W]
                Pout = Sn_rated;   
            
                % Copper losses [W]
                Pcu = (Ip^2 * R1 + Is^2 * R2); 

                % Input power [W]
                Pin = Pout + Ploss_steel + Pcu;  

                % Efficiency
                eta = Pout / Pin; 

                if eta >= Eff_required(i) && Cost_tot < min_cost
                    min_cost = Cost_tot;
                    best_B_core = B_core;
                    best_J_winding = J_winding;
                    best_K_insulation = K_insulation;
                    best_eta = eta;
                end
            end
        end
    end
end

% Print the best results
fprintf('Best cost: %f\n', min_cost);
fprintf('With B_core = %.2f T, J_winding = %.2f A/m^2, K_insulation = %.2f V\n', best_B_core, best_J_winding, best_K_insulation);