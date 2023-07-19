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

M = length(epochs);     % n° of observations = 280

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

% place the needed results from the integraton into X = (x, y, u, v, GM, Cd)
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
ylabel('$y_{\rho}$ $[km]$','interpreter','latex','FontSize',14)
title('Range Pre-Fit Residuals')

subplot(2,1,2)
plot(epochs,res_rhodot, 'x', 'color', '#009fd4')
xlabel('t\,[s]','Interpreter','latex','FontSize',12)
ylabel('$y_{\dot{\rho}}$ $[km/s]$','interpreter','latex','FontSize',14)
title('Range Rate Pre-Fit Residuals')

%% Building the Filter

% initializing the needed variables
max_iterations = 10;

rho_comp = zeros(M, 1);
rhodot_comp = zeros(M, 1);

X0_it = X0;
X0_it_save = zeros(max_iterations, N);
x_hat = 0;
x_bar = zeros(N, 1);

% list of apriori information
sigma_x = 0.5;          % km
sigma_y = 0.5;          % km
sigma_u = 7e-5;         % km/s
sigma_v = 7e-5;         % km/s
sigma_GM = 0.5;         % km^3/s^2
sigma_Cd = 0.15;

sigma_rho = 1*10^-3;         % km
sigma_rhodot = 50*10^-9;      % km/s

W_apr = eye(length(X0))*diag([1/sigma_x^2; 1/sigma_y^2; 1/sigma_u^2; 1/sigma_v^2; 1/sigma_GM^2; 1/sigma_Cd^2]);

% building of the W_obs matrix which contains the weights of each
% observables alternately repeated on the diagonal
W_obs = zeros(2*M);  

for i = 1 : M
    W_obs(2*i-1,2*i-1) = 1/(sigma_rho)^2; 
    W_obs(2*i,2*i) = 1/(sigma_rhodot)^2; 
end

% Start of the Iterative Process
for counter = 1 : max_iterations

    if counter >= 2
       for i = 1 : M
           W_obs(2*i-1,2*i-1) = 1/(RMS_rho)^2;
           W_obs(2*i,2*i) = 1/(RMS_rhodot)^2;
       end
       sigma_rho = RMS_rho;
       sigma_rhodot = RMS_rhodot;
    end
    
    PHI = eye(N);                   % at each cycle we reinitialize PHI
    phi = reshape(PHI, N*N, 1);  

    % now, we'll integrate the STM by giving as input the initial state
    % X0_it and phi as part of a single matrix [X0_it ; phi]

    [t, w] = ode113('Model_Transition', tspan, [X0_it ; phi], options);
    
    % then we must build H for each measurement 2:M+1
    H = [];
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
        H_t = H_tilde(X,Xs);
        
        H_j = H_t * PHI;
        
        H = [H ; H_j];      % we dock the support H_j underneath H
        
        % finally, we find the computed observables
        rho_comp(i-1) = sqrt((x-xs)^2 + (y-ys)^2);
        rhodot_comp(i-1) = ((x-xs)*u + (y-ys)*v)/rho_comp(i-1);
        
    end

    % now we evaluate the residuals
    y_rho = rho_obs - rho_comp;
    y_rhodot = rhodot_obs - rhodot_comp;

    % we define the epsilon vector
    eps = zeros(2*M, 1);
    
    % and we fill it with the values in order
    for i=1:length(epochs)
        eps(2*i-1) = y_rho(i);
        eps(2*i) = y_rhodot(i);
    end
    
    % WLS Correction
    P_cov = inv(H' * W_obs * H + W_apr);
    x_hat = (H' * W_obs * H + W_apr)\(H' * W_obs * eps + W_apr * x_bar);
    
    % State Update
    X0_it = X0_it +x_hat;
    X0_it_save(counter,:) = X0_it';
    x_bar = x_bar - x_hat;
    
    % Introduction of Parameters to check convergence
    SOS_rho = y_rho'*(1/sigma_rho)^2*y_rho;
    SOS_rhodot = y_rhodot'*(1/(sigma_rhodot)^2)*y_rhodot;

    MEAN_rho = mean(y_rho);
    MEAN_rhodot = mean(y_rhodot);
    
    RMS_rho = sqrt((y_rho'*y_rho)/length(y_rho));
    RMS_rhodot = sqrt((y_rhodot'*y_rhodot)/length(y_rhodot));
    
    % we need SOS similar to M = lenght(epochs)
    % and RMS similar to sigma_obs
    
    % Plotting the results with convergence parameters    
    figure(counter + 1)
    
    subplot(2,1,1)
    plot(epochs,y_rho,'x','color','#eb8900')
    xlabel('Seconds after $t_0= 0$','Interpreter','latex','FontSize',12)
    ylabel('$y_{\rho}$ [km]','interpreter','latex','FontSize',14)
    title(sprintf('Range Residuals at iteration n° %i \n SOS: %1.3e \t MEAN: %1.3e \t RMS: %1.3e',[counter,SOS_rho,MEAN_rho,RMS_rho]))
    
    
    subplot(2,1,2)
    plot(epochs,y_rhodot,'x','color','#009fd4')
    xlabel('Seconds after $t_0= 0$','Interpreter','latex','FontSize',12)
    ylabel('$y_{\dot{\rho}}$ [km/s]','interpreter','latex','FontSize',14)
    title(sprintf('Range Rate Residuals at iteration n° %i \n SOS: %1.3e \t MEAN: %1.3e \t RMS: %1.3e',[counter,SOS_rhodot,MEAN_rhodot,RMS_rhodot]))
    
end

counter = counter + 1;      % useful to continue naming figures

format longg
uncertainties = sqrt(diag(P_cov));
X0_f = X0_it_save(end, :);
final_table = [X0_f; uncertainties']

%% Investigation for most valuable parameter

partials_rho = zeros(M, N);
partials_rhodot = zeros(M, N);

for i=2:length(t)    
    H_til = H_tilde(X, Xs);
    H_t = H_til * PHI;    
    partials_rho(i-1, :) = H_t(1, :);    
    partials_rhodot(i-1, :) = H_t(2, :);    
end

% finally, we build the Information matrix for each observable
rho_info = partials_rho' * 1/sigma_rho^2 * partials_rho;
rhodot_info = partials_rhodot' * 1/sigma_rhodot^2 * partials_rhodot;

% the values will be on the diagonal for both range and rangerate
fprintf('Influence of Range on x_0: %3.3e\n',rho_info(1,1));
fprintf('Influence of Range on y_0: %3.3e\n',rho_info(2,2));
fprintf('Influence of Range on u_0: %3.3e\n',rho_info(3,3));
fprintf('Influence of Range on v_0: %3.3e\n',rho_info(4,4));
fprintf('Influence of Range on GM_0: %3.3e\n',rho_info(5,5));
fprintf('Influence of Range on Cd_0: %3.3e\n\n',rho_info(6,6));

fprintf('Influence of Range Rate on x_0: %3.3e\n',rhodot_info(1,1));
fprintf('Influence of Range Rate on y_0: %3.3e\n',rhodot_info(2,2));
fprintf('Influence of Range Rate on u_0: %3.3e\n',rhodot_info(3,3));
fprintf('Influence of Range Rate on v_0: %3.3e\n',rhodot_info(4,4));
fprintf('Influence of Range Rate on GM_0: %3.3e\n',rhodot_info(5,5));
fprintf('Influence of Range Rate on Cd_0: %3.3e\n\n',rhodot_info(6,6));

% for G = 1:6
%     ggg(G) = rho_info(G, G) > rhodot_info(G, G);
% end

fprintf(['We can see that the Range Rate influence on the initial state vector is larger for every parameter, therefore Range Rate\n' ...
    'is the most valuable observable.\n'])

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
    z_apo_min = 400;        % km

    if z_apo(j) < z_apo_min
        time = t24h(j);
        mnvr_time = seconds(time);
        mnvr_time.Format='hh:mm:ss';
        break
    end
end

fprintf('\nThe time at which we must perform the maneuver is:\n')
disp(mnvr_time)