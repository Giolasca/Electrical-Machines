function [A, d] = geom_core(Vp, N1, f, B_core)  
    
    % Cross-sectional area of the core based on flux density [m^2]
    A = Vp / (sqrt(2) * pi * N1 * f * B_core);
    
    % Diameter of the core from the cross-sectional area [m]
    d = sqrt(4 * A / pi);
    
end


