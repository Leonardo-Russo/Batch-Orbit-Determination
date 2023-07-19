% DEFINIZIONE DELLE EQUAZIONI DINAMICHE E LA DERIVATA DELLA STATE TRANSITION MATRIX

function dvect=Model_and_STM(t,vect)

% Parametri di Venere e dello spacecraft

RV=6052;
rho0=0.164;
z0=200;
QS=14.7;
area=40*10^(-6);
mass=2000;

dvect=zeros(42,1);

% Integrazione del vettore di stato

X=vect(1:6);
x=vect(1);
y=vect(2);
dx=vect(3);
dy=vect(4);
GM=vect(5);
Cd=vect(6);

r=sqrt(x^2+y^2);

% Componenti di Xdot

dvect(1)=dx;
dvect(2)=dy;
dvect(3)=(-GM/r^3)*x-1/2*Cd*rho0*(exp((RV+z0-r)/QS))*(area/mass)*sqrt((dx)^2+(dy)^2)*dx;
dvect(4)=(-GM/r^3)*y-1/2*Cd*rho0*(exp((RV+z0-r)/QS))*(area/mass)*sqrt((dx)^2+(dy)^2)*dy;
dvect(5)=0;
dvect(6)=0;

% Integrazione della state transition matrix

phi=vect(7:42);
PHI=reshape(phi,6,6);
A=A_matrix(X);
dphi=A*PHI;
dphi=reshape(dphi,1,36);

%Ccomponenti di phi_dot

dvect(7:42)=dphi;

end 