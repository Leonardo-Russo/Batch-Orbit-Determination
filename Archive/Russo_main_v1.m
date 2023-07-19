%% Homework 1 - Leonardo Russo 2015563
%% Workspace Initialization

clear all
close all
clc

%% Loading Observables and Defining Initial Quantities

% observables reading
obsData = load("observables.txt");
epochs = obsData(:, 1);         % s
rho_obs = obsData(:, 2);        % km
rhodot_obs = obsData(:, 3);     % km/s

M = length(epochs);     % n° of observations = 208

% initial condition
x0 = -0.8;          % km
y0 = 6419.4;        % km
u0 = -7.11389;      % km/s
v0 = -0.24912;      % km/s
GM = 324860.3;      % km^3/s^2
Cd = 2.2;

X0 = [x0; y0; u0; v0; GM; Cd];
X0_save = X0;
N = length(X0);     % n° of state variables = 6

Re = 6378;              % km
Rv = 6052;              % km
xe = -38.2e6;           % km 
xs = xe + Re*cosd(30);  % km
ys = Re*sind(30);       % km

Xs = [xs, ys];

% state vector X initialization
X = zeros(6, 1);

%% First Integration of the Trajectory

tspan = [0; epochs];    % we add the 0 since we integrate from t0 = 0

% defining the options for the integration
Tol0 = 1e-13;
Tol1 = 1e-13;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1);

[t, sol] = ode113('DynamicalModel2D', tspan, X0, options);

%_Plot of the Trajectory
%traj_plot(t, sol)

% place the integrated results into X

x = sol(2:end, 1);      % remember that length(x) == M
y = sol(2:end, 2);
u = sol(2:end, 3);
v = sol(2:end, 4);

% Pre-fit initialization and computation of computed observables

rho_comp_pre = zeros(M, 1);
rhodot_comp_pre = zeros(M, 1);

for i = 1 : M
    rho_comp_pre(i) = sqrt((x(i)-xs)^2+(y(i)-ys)^2);
    rhodot_comp_pre(i) = ((x(i)-xs)*u(i)+(y(i)-ys)*v(i))/rho_comp_pre(i);
end

res_rho = rho_obs - rho_comp_pre;
res_rhodot = rhodot_obs - rhodot_comp_pre;


% Plot of the pre-fit residuals
figure(1)

subplot(2,1,1)
plot(epochs,res_rho, 'x', 'color', '#eb8900')
xlabel('t\,[s]', 'Interpreter','latex','FontSize',12)
ylabel('$\epsilon_{\rho}$ $[km]$','interpreter','latex','FontSize',14)
title('Range Pre-Fit Residuals')

subplot(2,1,2)
plot(epochs,res_rhodot, 'x', 'color', '#009fd4')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$\epsilon_{\dot{\rho}}$ $[km/s]$','interpreter','latex','FontSize',14)
title('Range Rate Pre-Fit Residuals')

%% Building the Filter

% initializing the needed variables
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

sigma_range = 1*10^-3;         % km
sigma_rangerate = 50*10^9;    % km/s

% initial covariance matrix definition
P0bar = eye(N) * diag([sigma_x^2, sigma_y^2, sigma_u^2, sigma_v^2, sigma_GM^2, sigma_Cd^2]);

% building of the R matrix for observables precision
R = zeros(2*M);     
% we'll be placing the values alternately on the diagonal
for i = 1 : M
    R(2*i-1,2*i-1) = (sigma_range)^2; 
    R(2*i,2*i) = (sigma_rangerate)^2; 
end

% computation of the inverse of R
Rinv = inv(R);

% Start of the Iterative Process
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
        
        PHI = reshape(phi, N, N);   % reshapes phi into STM

        % we feed the values of H_t to build matrix H
        H_t = H_tilde(X, Xs);
        H(2*i-1:2*i, :) = H_t * PHI;

        % finally, we find the computed observables
        rho_comp(i-1) = sqrt((x-xs)^2 + (y-ys)^2);
        rhodot_comp(i-1) = ((x-xs)*u + (y-ys)*v)/rho_comp(i-1);

    end

    H = H(3:end, :);    % we remove the 2 initial rows in excess

    % now we evaluate the residuals
    y_rho = rho_obs - rho_comp;
    y_rhodot = rhodot_obs - rhodot_comp;
    
    % we define the epsilon vector
    y = zeros(2*M, 1);
    
    % and we fill it with the values in order
    for i = 1:length(epochs)
        y(2*i-1) = y_rho(i);
        y(2*i) = y_rhodot(i);
    end
    

    % Minimum Variance Estimate
    x_bar = zeros(N, 1);
    Pcov = inv(H' * Rinv * H + inv(P0bar));
    x_hat = (H' * Rinv * H + inv(P0bar))\(H' * Rinv * y +(inv(P0bar))*x_bar);
    
    % State Update
    x_bar = x_bar - x_hat;
    X0_it = X0_it + x_hat;
    
    % saving the current iteration value
    X0_it_save(counter,:) = X0_it;
    
    % Introduce Parameters to check Convergence
    
    SOS_rho = y_rho'*(1/sigma_range)^2*y_rho;
    SOS_rhodot = y_rhodot'*(1/(sigma_rangerate)^2)*y_rhodot;

    MEAN_rho = mean(y_rho);
    MEAN_rhodot = mean(y_rhodot);
    
    RMS_rho = sqrt((y_rho'*y_rho)/length(y_rho));
    RMS_rhodot = sqrt((y_rhodot'*y_rhodot)/length(y_rhodot));
    
    % we need SOS similar to M = lenght(epochs)
    % and RMS similar to sigma_obs
    
    %Grafici di controllo per la convergenza
    
    figure(counter + 1)
    
    subplot(2,1,1)
    plot(epochs,y_rho,'x','color','#eb8900')
    xlabel('Seconds after $t_0= 0$','Interpreter','latex','FontSize',12)
    ylabel('$\epsilon_{\rho}$ [km]','interpreter','latex','FontSize',14)
    title(sprintf('Range Residuals at iteration n° %i \n SOS: %1.3e \t MEAN: %1.3e \t RMS: %1.3e',[counter,SOS_rho,MEAN_rho,RMS_rho]))
    
    
    subplot(2,1,2)
    plot(epochs,y_rhodot,'x','color','#009fd4')
    xlabel('Seconds after $t_0= 0$','Interpreter','latex','FontSize',12)
    ylabel('$\epsilon_{\dot{\rho}}$ [km/s]','interpreter','latex','FontSize',14)
    title(sprintf('Range Rate Residuals at iteration n° %i \n SOS: %1.3e \t MEAN: %1.3e \t RMS: %1.3e',[counter,SOS_rhodot,MEAN_rhodot,RMS_rhodot]))
    
end

counter = counter + 1;

%% Investigation for most valuable parameter

partials_rho = zeros(M, N);
partials_rhodot = zeros(M, N);

for i=2:length(t)
    
    %PHI = reshape(s(i,6:30),5,5);
    
    H_til = H_tilde(X, Xs);
    
    H_t = H_til*PHI;
    
    partials_rho(i-1,:) = H_t(1,:);
    
    partials_rhodot(i-1,:) = H_t(2,:);
    
end

info_r = partials_rho'*partials_rho;

info_rr = partials_rhodot'*partials_rhodot;

fprintf('Influence of rho over x_0: %3.3e\n',info_r(1,1));
fprintf('Influence of rho over y_0: %3.3e\n',info_r(2,2));
fprintf('Influence of rho over u_0: %3.3e\n',info_r(3,3));
fprintf('Influence of rho over v_0: %3.3e\n',info_r(4,4));
fprintf('Influence of rho over GM_0: %3.3e\n\n',info_r(5,5));

fprintf('Influence of rho over x_0: %3.3e\n',info_rr(1,1));
fprintf('Influence of rho over y_0: %3.3e\n',info_rr(2,2));
fprintf('Influence of rho over u_0: %3.3e\n',info_rr(3,3));
fprintf('Influence of rho over v_0: %3.3e\n',info_rr(4,4));
fprintf('Influence of rho over GM_0: %3.3e\n',info_rr(5,5));


%% Evolution of the Semimajor axis and Eccentricity in 24 hours

% definition of 24h time span
t24h = 1:50:86400;
M24 = length(t24h);

% 24h trajectory integration using Dynamical Model
[t, g] = ode113('DynamicalModel2D', t24h, X0, options);

Rx = g(:, 1);       % km
Ry = g(:, 2);       % km
Vx = g(:, 3);       % km/s
Vy = g(:, 4);       % km/s
GMr = g(:, 5);      % km^3/s^2
Cdr = g(:, 6);

R_vect = [Rx'; Ry' ; zeros(1, M24)]';
V_vect = [Vx'; Vy'; zeros(1, M24)]';
R = zeros(M24, 1);
h_vect = zeros(M24, 3);
h = zeros(M24, 1);
e_vect = zeros(M24, 3);
e = zeros(M24, 1);
a = zeros(M24, 1);

% Computation of required parameters from ode results

for j = 1 : M24
    h_vect(j, :) = cross(R_vect(j, :), V_vect(j, :));
    R(j) = norm(R_vect(j, :));
    h(j) = norm(h_vect(j, :));

    e_vect(j, :) = ((cross(V_vect(j, :), h_vect(j, :)))/GM - (R_vect(j, :)/R(j)));
    e(j) = norm(e_vect(j, :));

    a(j) = ((h(j))^2/GM)*(1/(1-(e(j))^2));
end

z = R - Rv;

% Final Plots

figure(counter+1)
counter=counter+1;

plot(t24h, a, '-.', 'color', '#bd0f4f','LineWidth',1.5)
xlabel('Seconds after $t_0= 0$','Interpreter','latex','FontSize',12)
ylabel('Major Semiaxis $a$', 'Interpreter','latex','FontSize',12)
title('Semimajor Axis Evolution in 24h')

figure(counter+1)
counter=counter+1;
plot(t24h, e, '-.', 'color', '#4cad6e', 'LineWidth', 1.5)
xlabel('Seconds after $t_0= 0$','Interpreter','latex','FontSize',12)
ylabel('Eccentricity $e$', 'Interpreter','latex','FontSize',12)
title('Eccentricity Evolution in 24h')

%% Evolution of Altitude in 24 hours

figure(counter+1)

plot(t24h, z, '.', 'color', '#3f44e0', 'LineWidth',1.5)
xlabel('Seconds after $t_0= 0$','Interpreter','latex','FontSize',12)
ylabel('Altitude $z$', 'Interpreter','latex','FontSize',12)
title('Altitude Evolution in 24h')

%% Time at which we must perform the first maneuver

for j = 1 : M24
    R_apo(j) = ((h(j))^2/GM) * (1/(1-e(j)));
    z_apo(j) = R_apo(j) - Rv;

    if z_apo(j)<400
        time=t24h(j);
        time400km=seconds(time);
        time400km.Format='hh:mm:ss';
        break
    end
end

time400km