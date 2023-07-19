function Htilde = H_tilde(X,Xs)

    x = X(1);
    y = X(2);
    u = X(3);
    v = X(4);
    
    xs = Xs(1);
    ys = Xs(2);
    
    rho = sqrt((x-xs)^2+(y-ys)^2);
    
    Htilde = zeros(2,6);
    
    
    % H tilde matrix value assignment

    % Row 1
    Htilde(1,1) = (x-xs)/rho;
    Htilde(1,2) = (y-ys)/rho;

    % Row 2    
    Htilde(2,1) = u/rho - ((x-xs) * (u*(x-xs) + v*(y-ys)))/rho^3;   
    Htilde(2,2) = v/rho - ((y-ys) * (u*(x-xs) + v*(y-ys)))/rho^3;    
    Htilde(2,3) = (x - xs)/rho;    
    Htilde(2,4) = (y - ys)/rho;

end