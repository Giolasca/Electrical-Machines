%% Transformer Design Parameters Initialization
V_primary_nom = 10000; % Nominal primary voltage in Volts
V_secondary_nom = 400; % Nominal secondary voltage in Volts
S_rated = 1100e3; % Rated apparent power in VA for the whole transformer
transformer_ratio = V_primary_nom / V_secondary_nom; % Transformer ratio
turns_ratio = transformer_ratio; % Turns ratio, assuming ideal transformer for initial calculation

%% Design Variables Ranges
B_core_range = 0.5:0.25:2; % Flux density in Tesla
J_winding_range = 0.25:0.25:2; % Current density in A/mm^2 (adjust max value as needed)
K_insulation_range = 0.25:0.5:3; % Volts per turn

%% Other Constants
rho_Cu = 1.68e-8; % Resistivity of copper in ohm meter
frequency = 50; % Frequency in Hz

%% Iteration Setup
best_efficiency = 0; % Placeholder for the highest efficiency found
best_design = struct('B_core', 0, 'J_winding', 0, 'K_insulation', 0, 'efficiency', 0);

%% Iterative Design Process
for B_core = B_core_range
    for J_winding = J_winding_range
        for K_insulation = K_insulation_range
            % Calculate transformer parameters for current design variables
            
            N_primary = round(V_primary_nom / K_insulation);
            N_secondary = round(N_primary / turns_ratio);
            
            A_core = (V_secondary_nom * N_secondary) / (4.44 * frequency * B_core);
            
            I_primary_rated = S_rated / (sqrt(3) * V_primary_nom);
            I_secondary_rated = S_rated / (sqrt(3) * V_secondary_nom);
            
            A_wire_primary = I_primary_rated / J_winding;
            A_wire_secondary = I_secondary_rated / J_winding;
            
            L_primary = 2 * N_primary * sqrt(A_core); 
            L_secondary = 2 * N_secondary * sqrt(A_core);
            R_primary = rho_Cu * L_primary / A_wire_primary;
            R_secondary = rho_Cu * L_secondary / A_wire_secondary;
            
            % Placeholder efficiency calculation for illustration
            efficiency = 1 - (R_primary + R_secondary); 
            
            % Update best design if current design has better efficiency
            if efficiency > best_efficiency
                best_efficiency = efficiency;
                best_design.B_core = B_core;
                best_design.J_winding = J_winding;
                best_design.K_insulation = K_insulation;
                best_design.efficiency = efficiency;
            end
        end
    end
end

%% Display Best Design Parameters
fprintf('Best Design Parameters:\n');
fprintf('B_core: %.2f T\n', best_design.B_core);
fprintf('J_winding: %.2f A/mm^2\n', best_design.J_winding);
fprintf('K_insulation: %.2f V/turn\n', best_design.K_insulation);
fprintf('Efficiency: %.4f\n', best_design.efficiency);
%%
% Sensitivity Analysis Settings
sensitivity_range = 0.05; % 5% variation for sensitivity analysis

% Perform sensitivity analysis on the best design
for delta = -sensitivity_range: sensitivity_range: sensitivity_range
    % Vary design variables by a small percentage (delta)
    B_core_sensitivity = best_design.B_core * (1 + delta);
    J_winding_sensitivity = best_design.J_winding * (1 + delta);
    K_insulation_sensitivity = best_design.K_insulation * (1 + delta);
    
    % Recalculate the transformer cost based on new variables
    % Note: You would need a function or formula to calculate cost
    cost = calculate_transformer_cost(B_core_sensitivity, J_winding_sensitivity, K_insulation_sensitivity);
    
    % If the new cost is lower, update the design variables
    if cost < best_design.cost
        best_design.B_core = B_core_sensitivity;
        best_design.J_winding = J_winding_sensitivity;
        best_design.K_insulation = K_insulation_sensitivity;
        best_design.cost = cost;
    end
end
%%
% Update core and ECD parameters for three-phase design
best_design_3ph = best_design; % Start with single-phase design

% Update the parameters that change in three-phase configuration
% Note: You will need to define the relationships between single-phase and three-phase parameters
best_design_3ph.A_core_3ph = update_core_area_3ph(best_design.A_core);
best_design_3ph.R_primary_3ph = update_winding_resistance_3ph(best_design.R_primary);
best_design_3ph.R_secondary_3ph = update_winding_resistance_3ph(best_design.R_secondary);
% ... and other parameters as necessary

% Display updated design parameters for three-phase transformer
fprintf('Three-Phase Design Parameters:\n');
fprintf('Core Area: %.2f mm^2\n', best_design_3ph.A_core_3ph);
fprintf('Primary Resistance: %.4f Ohms\n', best_design_3ph.R_primary_3ph);
fprintf('Secondary Resistance: %.4f Ohms\n', best_design_3ph.R_secondary_3ph);
% ... and other parameters as necessary
