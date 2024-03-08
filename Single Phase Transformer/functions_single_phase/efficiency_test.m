function efficiency = efficiency_test(V_p, V_s, S_n, R_1, X_1, R_2, X_2, P_steel, flag)
    
    % Define the vector of power factors
    PF = linspace(0.1,0.9,9);
    S = zeros(length(PF),1);
    
    % Calculate the new powers
    for i = 1:length(PF)
        if flag
            S(i) = S_n * (PF(i) + 1i*sqrt(1 - PF(i)^2));
        else
            S(i) = S_n * (PF(i) - 1i*sqrt(1 - PF(i)^2));
        end
    end

    I1 = S / V_p;
    I2 = S / V_s;
    Z1 = R_1 + 1i*X_1;
    Z2 = R_2 + 1i*X_2;
    
    Pout = S;
    Pcu = Z1 * I1.^2 + Z2 * I2.^2;
    Pin = Pcu + P_steel + Pout;
    
    % Efficiency
    efficiency = abs(Pout) ./ abs(Pin);
end

