function conf = save_conf(params, geometry, electrical, costs)
    conf = struct(...
        'Parameters', struct(...
            'B_core', params(1), ...
            'J_winding', params(2), ...
            'K_insulation', params(3), ...
            'h_windings', params(4)...
        ), ...
        'Geometry', struct(...
            'Section_wire_primary', geometry(1), ...
            'Lenght_wire_primary', geometry(2), ...
            'Volume_wire_primary', geometry(3), ...
            'Section_wire_secondary', geometry(4), ...
            'Lenght_wire_secondary', geometry(5), ...
            'Volume_wire_secondary', geometry(6), ...
            'Section_core', geometry(7), ...
            'heigth_core', geometry(8), ...
            'width_core', geometry(9), ...
            'Volume_core', geometry(10)...
        ), ...
        'Electrical', struct(...
            'R1p', electrical(1), 'X1p', electrical(2), 'Z1p', electrical(3), ...
            'R2p', electrical(4), 'X2p', electrical(5), 'Z2p', electrical(6), ...
            'R0', electrical(7), 'Xm', electrical(8), 'Zm', electrical(9)...
        ), ...
        'Costs', struct(...
            'Cost_steel', costs(1), ...
            'Cost_copper', costs(2), ...
            'Cost_oil', costs(3), ...
            'Cost_tank', costs(4), ...
            'Total_Cost', costs(5)...
        )...
    );
end
