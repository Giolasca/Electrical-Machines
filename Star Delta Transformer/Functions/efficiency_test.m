function [performance_table] = efficiency_test(Vpp, S_n, Z_1p, Z_2p, Zm, P_steel)
   
    % Define the vector of power factors
    PF_descending = linspace(0.9, 0, 10);
    PF_ascending = linspace(0.1, 0.9, 9);
    PF = [PF_descending, PF_ascending(1:end)];  

    S = zeros(length(PF), 1);
    
    % Calculate the new powers
    for i = 1:length(PF)
        if (i >= 0)
            S(i) = S_n * (PF(i) + 1i * sqrt(1 - PF(i)^2));
        else
            S(i) = S_n * (PF(i) - 1i * sqrt(1 - PF(i)^2));
        end
    end
    
    % P_out [W]
    Pout = real(S);

    % Close the circuit with load
    RLp = Vpp^2 ./ Pout;
    XLp = Vpp^2 ./ imag(S);

    % Total load impedence
    ZLp = RLp + 1i*XLp; 

    % Primary current from equivalent circuit
    Zeq_power = Z_1p + ((Z_2p + ZLp) .* Zm) ./ (Zm + Z_2p + ZLp);

    % Equivalent current referred to primary [A]
    I1_power = Vpp ./ Zeq_power;

    % Power [W]
    S1_power = Vpp .* conj(I1_power);

    % Active input power related to copper [W]
    Pcu = real(S1_power);

    % Input power [W]
    Pin = Pout + P_steel + Pcu;

    % Efficiency
    Efficiency = Pout ./ Pin;
        
    % Voltage regulation
    I_rated = S ./(Vpp);
    V1_power = Vpp + (Z_1p + Z_2p) .* conj(I_rated);
    Voltage_regulator = (abs(V1_power) - Vpp)*100/Vpp;
    
    % Create the table
    performance_table = table(PF', Efficiency, Voltage_regulator, S, I_rated, V1_power, ...
                             'VariableNames', {'Power_factor', 'Efficiency', 'Voltage_regulator', 'S', 'I_rated', 'V1_power'});

end