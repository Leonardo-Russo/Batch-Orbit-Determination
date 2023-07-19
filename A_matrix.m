function [A] = A_matrix(X)
% Input: X is the state vector
% Purpose: this function computes the Amatrix for the given state and model
    
    % Initialization and variables assignment
    A = zeros(6,6);

    x = X(1);
    y = X(2);
    u = X(3);
    v = X(4);
    GM = X(5);
    Cd = X(6);

    r = sqrt(x^2+y^2);
    V = sqrt(u^2+v^2);

    S = 40e-6;              % km^2
    m = 2000;               % kg
    H = 14.7;               % km

    % Rows 1-2
    A(1, 3) = 1;
    A(2, 4) = 1;

    % Row 3
    A(3, 1) = -GM/r^3 + (3*GM*x^2)/r^5 + 0.5*Cd*S/m*V*u*x/(H*r)*elvtn_density(r);
    A(3, 2) = 3*GM*x*y/r^5 + 0.5*Cd*S/m*V*u*y/(H*r)*elvtn_density(r);
    A(3, 3) = -0.5*Cd*S/m*elvtn_density(r)*(V+u^2/V);
    A(3, 4) = -0.5*Cd*S/m*u*elvtn_density(r)*v/V;
    A(3, 5) = -x/r^3;
    A(3, 6) = -0.5*elvtn_density(r)*S/m*V*u;

    % Row 4
    A(4, 1) = 3*GM*x*y/r^5 + 0.5*Cd*S/m*V*v*x/(H*r)*elvtn_density(r);
    A(4, 2) = -GM/r^3 + 3*GM*y^2/r^5 + 0.5*Cd*S/m*V*v*y/(H*r)*elvtn_density(r);
    A(4, 3) = -0.5*Cd*S/m*v*u/V*elvtn_density(r);
    A(4, 4) = -0.5*Cd*S/m*elvtn_density(r)*(V+v^2/V);
    A(4, 5) = -y/r^3;
    A(4, 6) = -0.5*elvtn_density(r)*S/m*V*v;

    % Rows 5-6
    % no need to define them since they would be zero anyway

end