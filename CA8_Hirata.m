%% Coding Assignment #8, MCDB 108C, Spring 2024 
% Description: OSCILLATORY SYSTEM
% Anthony Sacco

clearvars
%% Parameters
% Define parameters %INSERT 1
A = 0.022;  
B = 0.3; 
C = 0.031; 
D = 0.028;
E = 0.5; 
F = 20;
G = 0.3;

% Define time vector and steps %INSERT 2
t_min = 0; %start time in minute
t_max = 1000; %end time in minute
dt = 0.01; %time step
t = 0:dt:(1000-dt); %time vector
n_steps = length(t); %number of time steps

% Create & initialize vectors for the normal condition %INSERT 3
X = zeros(1, n_steps);
Y = zeros(1, n_steps);
Z = zeros(1, n_steps);

% Create & initialize vectors for Type I perturbation %INSERT 3
X1 = zeros(1, n_steps);
Y1 = zeros(1, n_steps);
Z1 = zeros(1, n_steps);

% Create & initialize vectors for Type II perturbation %INSERT 3
X2 = zeros(1, n_steps);
Y2 = zeros(1, n_steps);
Z2 = zeros(1, n_steps);

% Create perturbation Type I %READ
perturbation = zeros(n_steps,1);
perturbation(ceil(270/dt)) = 10; % Here ceil is used to unsure that the result of the division 270/dt is a whole number that can be used as an index


%% Part 1: Euler's integration of ODE system with and without perturbations
for i = 1:n_steps-1

    % Normal condition %INSERT 5
    dXdt = (B*Y(i)) - (A*X(i)*Z(i)) - (C*X(i));  % Hes1 protein
    dYdt = (E/(1+X(i)^2)) - (D*Y(i));  % Hes mRNA
    dZdt = (F/(1+X(i)^2)) - (A*X(i)*Z(i)) - (G*Z(i));  % Hes1-interacting factor

    X(i+1) = X(i) + dXdt * dt;
    Y(i+1) = Y(i) + dYdt * dt;
    Z(i+1) = Z(i) + dZdt * dt;

    % Perturbation I %INSERT 6
    if perturbation(i) ~= 0
        Y1(i) = Y1(i) + perturbation(i);
    end

    dXdt_perturb1 = (B*Y1(i)) - (A*X1(i)*Z1(i)) - (C*X1(i));  % Hes1 protein
    dYdt_perturb1 = (E/(1+X1(i)^2)) - (D*Y1(i));  % Hes mRNA
    dZdt_perturb1 = (F/(1+X1(i)^2)) - (A*X1(i)*Z1(i)) - (G*Z1(i));  % Hes1-interacting factor

    X1(i+1) = X1(i) + dXdt_perturb1 * dt;
    Y1(i+1) = Y1(i) + dYdt_perturb1 * dt;
    Z1(i+1) = Z1(i) + dZdt_perturb1 * dt;

    % Perturbation II %INSERT 7
    if i == ceil(300/dt)
        Y2(i) = 0; % Create perturbation at time t = 300 min
    elseif i == ceil(600/dt)
        Y2(i) = 0; % Create perturbation at time t = 600 min
    end

    dXdt_perturb2 = (B*Y2(i)) - (A*X2(i)*Z2(i)) - (C*X2(i));  % Hes1 protein
    dYdt_perturb2 = (E/(1+X2(i)^2)) - (D*Y2(i));  % Hes mRNA
    dZdt_perturb2 = (F/(1+X2(i)^2)) - (A*X2(i)*Z2(i)) - (G*Z2(i));  % Hes1-interacting factor
    
    X2(i+1) = X2(i) + dXdt_perturb2 * dt;
    Y2(i+1) = Y2(i) + dYdt_perturb2 * dt;
    Z2(i+1) = Z2(i) + dZdt_perturb2 * dt;

end 

%% Part 2: Plotting and analyzing the results
close all;

figure (1);
subplot (3,1,1)
plot(t,X, 'g'); hold on;
plot(t,Y,'r');
plot(t,Z,'k');
legend ('X (hes1 protein)', 'Y (hes1 mRNA)','Z (hes1 interacting factor')
title('Oscillations of hes1 mRNA, protein, and interaction factor')

subplot (3,1,2)
plot (t, X1, 'g'); hold on
plot (t, Y1, 'r');
plot (t, Z1, 'k');
ylabel('Expression level', 'FontSize',16)
title('Oscillations with perturbation (type I)')

subplot (3,1,3)
plot (t, X2, 'g'); hold on
plot (t, Y2, 'r');
plot (t, Z2, 'k');
xlabel('Time (min)','FontSize',16)
title('Oscillations with perturbation (type II)')

figure (2);
subplot (3,1,1)
plot(X,Y); 
title ('2D trajectory of system ')

subplot(3,1,2)
plot(X1, Y1)
ylim([0 4])
ylabel ('Y (hes1 mRNA)', 'FontSize', 16)
title ('2D trajectory of system with perturbation type I')
 
subplot(3,1,3)
plot(X2, Y2)
xlabel ('X (hes1 protein)', 'FontSize', 16)
title ('2D trajectory of system with perturbation type II')

figure (3);
subplot (3,1,1)
plot3(X,Y,Z); grid on;
title ('3D trajectory of system ')

subplot(3,1,2)
plot3(X1, Y1, Z1); grid on;
title ('3D trajectory of system with perturbation type I')
 
subplot(3,1,3)
plot3(X2, Y2, Z2); grid on;
xlabel ('X (hes1 protein)', 'FontSize', 16)
ylabel ('Y (hes1 mRNA)', 'FontSize', 16)
zlabel ('Z (hes1 interacting factor', 'FontSize', 16)
title ('3D trajectory of system with perturbation type II')

%% COMMENT: Comment your observation about the system under normal condition vs. with the 2 types of perturbations 

% For the normal condition the system shows synchronized oscillations for
% the Hes1 protein, the Hes1 mRNA, and the interacting factors. However,
% looking at the type 1 perturbation the system displays different behavior.
% Upon perturbation(increasing Hes1 mRNA Y1) at 270 minutes there were 
% significant spikes for the Hes1 protein(X1) and the interacting factor(Z1). This results in a
% disruption of the Notch signaling pathway and its negative feedback.
% Due to this increase in the Hes1 protein, the negative feedback is still
% allowed to act which is why you eventually see the system return back to
% the synchronized oscillations(equilibrium) see under the normal
% condition. The type 2 perturbation of the system also displayed different
% behavior compared to the normal condition and the type 1 perturbation.
% This is because Hes1 mRNA was effectively forced to zero at time 300 and
% 600 resulting in the decrease of Hes1 protein(X2), Hes1 mRNA(Y2), and the
% interacting factor(Z2). This once again results in the disruption of the
% Notch signaling pathway. Since the concentration of Hes1 protein
% decreases after each perturbation, there is not enough protein to engage
% the negative feedback loop which is why you see an increase in Hes1
% protein back to equilibrium after each perturbation.