function rho = elvtn_density(r)
% Input: distance from Venus center r in [km]

    % definition of rho and costants
    Rv = 6052;          % km
    z = r - Rv;         % km

    z0 = 200;           % km
    rho0 = 1.64e-1;     % kg/km^3
    H = 14.7;           % km

    % now we apply the formula

    rho = rho0 * exp(-(z-z0)/H);

end

    
    