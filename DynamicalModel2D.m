function ds = DynamicalModel2D(t, s)
% Input: s = (x, y, u, v, GM, Cd)
    
    % we initialize the vector and define quantities for Xdot
    ds = zeros(6, 1);

    x = s(1);
    y = s(2);
    u = s(3);
    v = s(4);
    GM = s(5);
    Cd = s(6);

    r = sqrt(x^2+y^2);
    V = sqrt(u^2+v^2);

    % constants - must be defined in the problem?
    A = 40e-6;              % km^2
    m = 2000;               % kg
    
    % value assignment
    ds(1) = u;
    ds(2) = v;
    ds(3) = -GM/r^3 * x - 0.5 * elvtn_density(r) * Cd * A/m * V * u;
    ds(4) = -GM/r^3 * y - 0.5 * elvtn_density(r) * Cd * A/m * V * v;
    
    % no need to define further because they would be zeros anyway

end

