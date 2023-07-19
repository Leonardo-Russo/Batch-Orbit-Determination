function ds = Model_Transition(t,s)

    ds = zeros(6+6*6,1);
    
    X = s(1:6);
    x = s(1);
    y = s(2);
    u = s(3);
    v = s(4);
    GM = s(5);
    Cd = s(6);
    
    r = sqrt(x^2 + y^2);
    V = sqrt(u^2+v^2);
    
    % problem constants
    A = 40e-6;              % km^2
    m = 2000;               % kg
    
    % State Integration
    ds(1) = u;
    ds(2) = v;
    ds(3) = -GM/r^3 * x - 0.5 * elvtn_density(r) * Cd * A/m * V * u;
    ds(4) = -GM/r^3 * y - 0.5 * elvtn_density(r) * Cd * A/m * V * v;
    % similarly to before, the last two just remain zeros

    % we can now compute matrix A with a previoulsly defined function
    A = A_matrix(X);
    
    % then we can integrate the State Transition Matrix PHI
    PHI = reshape(s((6+1):end),6,6);    % we reshape from N+1
    
    dPHI = A*PHI;
    
    ds((6+1):end) = reshape(dPHI,6*6,1);

end