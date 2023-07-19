% FUNZIONE CHE DEFINISCE LA MATRICE A (DINAMICA)

function [A]=MatriceA(X)

% Parametri di Venere e dello spacecraft

RV=6052;
rho0=0.164;
z0=200;
QS=14.7;
area=40*10^(-6);
mass=2000;

A = zeros (6,6);

x  = X(1);
y  = X(2);
u  = X(3);
v  = X(4);
GM = X(5);
Cd = X(6);

r=sqrt(x^2+y^2);
z=r-RV;
rhoz= rho0*exp( -(z-z0)/QS);

% Componenti della matrice A

A(1,3) = 1;
A(2,4) = 1;

A(3,1) = -(GM/r^3)+(3*GM*x^2/(r^5))   -   (1/2)*Cd*(area/mass)*sqrt(u^2+v^2)*u*(-rho0*x*exp((RV+z0-r)/QS)/(QS*r));
A(3,2) =  3*GM*x*y/(r^5)   -   (1/2)*Cd*(area/mass)*sqrt(u^2+v^2)*u*(-rho0*y*exp((RV+z0-r)/QS)/(QS*r) ) ;
A(3,3) = -(1/2)*rhoz*Cd*(area/mass)*sqrt(u^2+v^2)   -   (1/2)*rhoz*Cd*(area/mass)*u*(u/sqrt(u^2+v^2));
A(3,4) = -(1/2)*rhoz*Cd*(area/mass)*u*(v/sqrt(u^2+v^2));
A(3,5) = -x/r^3;
A(3,6) = -(1/2)*rhoz*(area/mass)*sqrt(u^2+v^2)*u;

A(4,1) =  3*GM*x*y/(r^5)   -   (1/2)*Cd*(area/mass)*sqrt(u^2+v^2)*v*(-rho0*x*exp((RV+z0-r)/QS)/(QS*r) ) ;
A(4,2) = -(GM/r^3)+(3*GM*y^2/(r^5))   -   (1/2)*Cd*(area/mass)*sqrt(u^2+v^2)*v*(-rho0*y*exp((RV+z0-r)/QS)/(QS*r) ) ;
A(4,3) = -(1/2)*rhoz*Cd*(area/mass)*u*(u/sqrt(u^2+v^2));
A(4,4) = -(1/2)*rhoz*Cd*(area/mass)*sqrt(u^2+v^2)   -   (1/2)*rhoz*Cd*(area/mass)*v*(v/sqrt(u^2+v^2));
A(4,5) = -y/r^3;
A(4,6) = -(1/2)*rhoz*(area/mass)*sqrt(u^2+v^2)*v;

end