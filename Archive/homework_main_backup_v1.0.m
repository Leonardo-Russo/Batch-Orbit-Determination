%% Initialization

clear all
close all
clc

%% Loading Observables and Defining Initial Values

obsData = load("observables.txt");
epochs = obsData(:, 1);         % s
rho_obs = obsData(:, 2);        % km
rhodot_obs = obsData(:, 3);     % km/s

M = length(epochs);     % n° of observations

% Initial Values

x0 = -0.8;          % km
y0 = 6419.4;        % km
u0 = -7.11389;      % km/s
v0 = -0.24912;      % km/s
GM = 324860.3;      % km^3/sec^2
Cd = 2.2;

X0 = [x0; y0; u0; v0; GM; Cd];
X0_save = X0;
N = length(X0);

Re = 6378;          % km
Rv = 6052;          % km
xe = -38.2e6;       % km 
xs = xe+Re;         % km
ys = Re*sind(30);   % km

Xs = [xs, ys];

% finally, we initialize the state vector X
X = zeros(6, 1);

%% Integrate the Trajectory

tspan = [0; epochs];    % we must add the 0 since we integrate from t0=0

% Defining the options for the integration
Tol0 = 1e-13;
Tol1 = 1e-13;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1);

% tspan_temp = [0, epochs(1)];    % is it useful?

[t, sol] = ode113('DynamicalModel2D', tspan, X0, options);

%##% additional info: plot of the trajectory
traj_plot(t, sol)

% considera -Rv

%% 
% we can place the values into X

x = sol(2:end, 1);
y = sol(2:end, 2);
u = sol(2:end, 3);
v = sol(2:end, 4);

%##% note that length(x) == M

% Computed observables initialization and computation

rho_comp_pre = zeros(M, 1);
rhodot_comp_pre = zeros(M, 1);

for i = 1 : M

    rho_comp_pre(i) = sqrt((x(i)-xs)^2+(y(i)-ys)^2);
    rhodot_comp_pre(i) = ((x(i)-xs)*u(i)+(y(i)-ys)*v(i))/rho_comp_pre(i);

end

%##% additional info: chech the size of residuals
%checkresiduals(rho_obs, rhodot_obs, rho_comp, rhodot_comp);

% GOT TO CHANGE HERE
res_rho = rho_obs - rho_comp_pre;
res_rhodot = rhodot_obs - rhodot_comp_pre;

figure(1)

subplot(2,1,1)
plot(epochs,res_rho,'b+')
title('Errore sul range')
xlabel('t [s]')
ylabel('$\epsilon_{\rho}$ $[km]$','interpreter','latex','FontSize',15)

subplot(2,1,2)
plot(epochs,res_rhodot,'b+')
title('Errore sul range rate')
xlabel('t [s]')
ylabel('$\epsilon_{\dot{\rho}}$ $[\frac{km}{s}]$','interpreter','latex','FontSize',15)

%%%% perchè mi dice che è iterativo?????????

%% Build the FUCKING Filter

% we initialize the needed variables
max_iterations = 6;

rho_comp = zeros(M, 1);
rhodot_comp = zeros(M, 1);
H = zeros(2*M, N);

X0_it = X0;
X0_it_save = zeros(max_iterations, N);

x_hat = 0;

% we define the apriori weight matrix --- IT MUST BE CHANGED!!!!!!!!
W_apr = eye(N)*diag([1/100;1/100;1/0.01;1/0.01;1/9]);

W_obs = zeros(2*M);

for i = 1 : M
    W_obs(2*i-1,2*i-1) = 1/(1e-4)^2; %Pesi per il range
    W_obs(2*i,2*i) = 1/(1e-3)^2; %Pesi per il range rate
end

% then we start the iteration
for counter = 1 : max_iterations

    PHI = eye(N);                   % at each cycle we initialize PHI
    phi = reshape(PHI, N*N, 1);     

    % no, we'll integrate the STM by giving as input the initial state
    % X0_it and phi as part of a single matrix [X0_it ; phi]

    [t, w] = ode113('Model_Transition', tspan, [X0_it ; phi], options);

    % then we must build H for each measurement 2:M

    for i = 2 : M
        
        X = w(i, 1:5);
        x = X(1);
        y = X(2);
        u = X(3);
        v = X(4);

        phi = w((N+1):(N*N+N));     % retrieves phi as a vector
        PHI = reshape(phi, N, N);   % reshapes phi into STM

        Htilde = H_tilde(X, Xs);
        H_t = Htilde * PHI;

        % now we put the values of H_it into the right spots of H
        H(2*(i-1)-1,:) = H_t(1,:);
        H(2*(i-1),:) = H_t(2,:);

        % finally, we find the computed observables
        rho_comp(i-1) = sqrt((x-xs)^2 + (y-ys)^2);
        rhodot_comp(i-1) = ((x-xs)*u + (y-ys)*v)/rho_comp(i-1);
        
        
        % orbiter code
        %rho_comp = sqrt((x-xs)^2+(y-ys)^2);
        %rhodot_comp(i-1) = 1/rho_comp * ((x-xs)*u*(y-ys)*v);

    end

    % now we evaluate the residuals
    eps_rho = rho_obs - rho_comp;
    eps_rhodot = rhodot_obs - rhodot_comp;
    
    % we define the epsilon matrix
    eps = zeros(2*M, 1);
    
    % and we place the correct values inside 
    for i=1:length(epochs)
        eps(2*i-1) = eps_rho(i);
        eps(2*i) = eps_rhodot(i);
    end
    
    %!!!!!% TO CHANGE

    %Propagazione delle informazioni a priori 
    x_hat_apr = X0 - X_0;
    
    x_bar = PHI*x_hat_apr;
    
    
    
    %MISSING

end





