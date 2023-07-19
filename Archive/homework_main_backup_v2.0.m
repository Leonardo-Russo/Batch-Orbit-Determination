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
GM = 324860.3;      % km^3/s^2
Cd = 2.2;

X0 = [x0; y0; u0; v0; GM; Cd];
X0_save = X0;
N = length(X0);

Re = 6378;              % km
Rv = 6052;              % km
xe = -38.2e6;           % km 
xs = xe + Re*cosd(30);    % km
ys = Re*sind(30);       % km

Xs = [xs, ys];

% finally, we initialize the state vector X
X = zeros(6, 1);

%% Integrate the Trajectory

tspan = [0; epochs];    % we must add the 0 since we integrate from t0=0

% Defining the options for the integration
Tol0 = 1e-13;
Tol1 = 1e-13;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1);

[t, sol] = ode113('DynamicalModel2D', tspan, X0, options);

%##% additional info: plot of the trajectory
%traj_plot(t, sol)

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
%checkresiduals(rho_obs, rhodot_obs, rho_comp_pre, rhodot_comp_pre);

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

%% Build the FUCKING Filter

% we initialize the needed variables
max_iterations = 6;

rho_comp = zeros(M, 1);
rhodot_comp = zeros(M, 1);

H = zeros(2*M, N);

X0_it = X0;
X0_it_save = zeros(max_iterations, N);

x_hat = 0;

% list of apriori information
sigma_x = 0.5;          % km
sigma_y = 0.5;          % km
sigma_u = 7e-5;         % km/s
sigma_v = 7e-5;         % km/s
sigma_GM = 0.5;         % km^3/s^2
sigma_Cd = 0.15;

sigma_range = 1e-3;         % km
sigma_rangerate = 50e-9;    % km/s

% building of the Weight matrix from apriori information
W_apr = eye(N) * diag([1/sigma_x^2, 1/sigma_y^2, 1/sigma_u^2, 1/sigma_v^2, 1/sigma_GM^2, 1/sigma_Cd^2]);

% building of the Weight matrix for observables
W_obs = zeros(2*M);     % we're placing the values on the diagonal alternated

for i = 1 : M

    W_obs(2*i-1,2*i-1) = 1/(sigma_range)^2; 
    W_obs(2*i,2*i) = 1/(sigma_rangerate)^2; 
    
end

% we start the Iterative Process
for counter = 1 : max_iterations

    PHI = eye(N);                   % at each cycle we reinitialize PHI
    phi = reshape(PHI, N*N, 1);     

    % now, we'll integrate the STM by giving as input the initial state
    % X0_it and phi as part of a single matrix [X0_it ; phi]

    [t, w] = ode113('Model_Transition', tspan, [X0_it ; phi], options);

    % then we must build H for each measurement 2:M+1

    for i = 2 : M+1
        
        X = w(i, 1:6);
        x = X(1);
        y = X(2);
        u = X(3);
        v = X(4);
        GM = X(5);
        Cd = X(6);

        phi = w(i, (N+1):end);     % retrieves phi as a vector
        % check that end == 36+6 = 42
        PHI = reshape(phi, N, N);   % reshapes phi into STM

        Htilde = H_tilde(X, Xs);
        H_t = Htilde * PHI;

        % we put the values of H_t into the right spots of H by row vectors
        H(2*(i-1)-1,:) = H_t(1,:);
        H(2*(i-1),:) = H_t(2,:);

        % finally, we find the computed observables
        rho_comp(i-1) = sqrt((x-xs)^2 + (y-ys)^2);
        rhodot_comp(i-1) = ((x-xs)*u + (y-ys)*v)/rho_comp(i-1);

    end

    % now we evaluate the residuals
    eps_rho = rho_obs - rho_comp;
    eps_rhodot = rhodot_obs - rhodot_comp;
    
    % we define the epsilon vector
    y = zeros(2*M, 1);
    
    % and we fill it with the values in order
    for i=1:length(epochs)
        y(2*i-1) = eps_rho(i);
        y(2*i) = eps_rhodot(i);
    end
    
    %!!!!!% TO CHANGE

    % Apriori Information Matrix Update
    x_hat_apr = X0_it - X0;
    x_bar = PHI * x_hat_apr;

    %Stima della correzione da apportare
    
    P_cov = inv(H' * W_obs * H + W_apr);
    
    x_hat = (H' * W_obs * H + W_apr)\(H' * W_obs * y + W_apr * x_bar);
    % this writing is faster and more accurate for Matlab
    
    % initial state correction
    X0_it = X0_it + x_hat;
    
    % save the current iteration value
    X0_it_save(counter,:) = X0_it';
    
    
    %Valutazione dei parametri per controllare la convergenza
    
    mean_range = mean(eps_rho);
    mean_rangerate = mean(eps_rhodot);
    
    RMS_range = sqrt((eps_rho'*eps_rho)/length(eps_rho));
    RMS_rangerate = sqrt((eps_rhodot'*eps_rhodot)/length(eps_rhodot));
    
    SOS_range = eps_rho'*(1/(sigma_range)^2)*eps_rho;
    SOS_rangerate = eps_rhodot'*(1/(sigma_rangerate)^2)*eps_rhodot;

    % we need SOS similar to M = lenght(epochs)
    
    
    
    %Grafici di controllo per la convergenza
    
    figure(counter + 1)
    
    subplot(2,1,1)
    plot(epochs,eps_rho,'b+')
    title(sprintf('Range residuals all''iterazione %i \n M: %1.4e \t RMS: %1.4e \t SOS: %1.4e',[counter,mean_range,RMS_range,SOS_range]))
    xlabel('Secondi dopo t_0=0')
    ylabel('$\epsilon_{\rho}$ [km]','interpreter','latex','FontSize',15)
    
    
    subplot(2,1,2)
    plot(epochs,eps_rhodot,'b+')
    title(sprintf('Range rate residuals all''iterazione %i \n M: %1.4e \t RMS: %1.4e \t SOS: %1.4e',[counter,mean_rangerate,RMS_rangerate,SOS_rangerate]))
    xlabel('Secondi dopo t_0=0')
    ylabel('$\epsilon_{\dot{\rho}}$ $\left[\frac{km}{s}\right]$','interpreter','latex','FontSize',15)
    

end





