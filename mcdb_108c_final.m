%% Final Coding Assignment
% Anthony Sacco
% This is the result of my own work along with referencing previous HW's as
% suggested by the instructions.

clearvars
%% Define the parameters

%Definte time vector and steps
K = 0.1;    % Threshold of sigmoid
n = 10;    % Hill coefficient
dt = 1 / 3600;   % Time step in seconds
t = 0:dt:(100 - dt); %Time vector in hours
num_steps = length(t); % Number of steps

% Initialize vectors with and without perturbation (x2 = perturbation,etc.)
x = zeros(1, num_steps);
y = zeros(1, num_steps);
z = zeros(1, num_steps);
x2 = zeros(1, num_steps);
y2 = zeros(1, num_steps);
z2 = zeros(1, num_steps);

% Initial conditions with and without perturbation
x(1) = 0;
y(1) = 0;
z(1) = 0;
x2(1) = 0;
y2(1) = 0;
z2(1) = 0;

%% Define the ODE system using Euler's method
%Integration without perturbation
for i = 1:num_steps-1
    dxdt = ((K^n) / (K^n + z(i)^n)) - x(i);
    dydt = x(i) - y(i);
    dzdt = y(i) - z(i);
    
    x(i+1) = x(i) + dxdt * dt;
    y(i+1) = y(i) + dydt * dt;
    z(i+1) = z(i) + dzdt * dt;
end

    %Integration with perturbation
for i = 1:num_steps-1
    if i == ceil(50/dt)
        x2(i) = 0.13;
        y2(i) = 0.13;
        z2(i) = 0.13;
    end

    dxdt_perturb = ((K^n) / (K^n + z2(i)^n)) - x2(i);
    dydt_perturb = x2(i) - y2(i);
    dzdt_perturb = y2(i) - z2(i);

    x2(i+1) = x2(i) + dxdt_perturb * dt;
    y2(i+1) = y2(i) + dydt_perturb * dt;
    z2(i+1) = z2(i) + dzdt_perturb * dt;
end

%% Plot the results
close all;

%Plot without perturbation
figure (1);
subplot(2,1,1);
plot(t, x, 'r', 'LineWidth', 2);
hold on;
plot(t, y, 'b', 'LineWidth', 2);
plot(t, z, 'g', 'LineWidth', 2);
xlabel('Time (Hours)');
ylabel('Concentration (arbitrary unit)');
legend('x (mRNA Concentration)','y (Protein Concentration)','z (Inhibiting Protein Concentration)');
title('Time course of x, y and z in absence of perturbation');
xlim([0 100]);  
ylim([0 0.7]);

%Plot with perturbation
subplot(2,1,2);
plot(t, x2, 'r', 'LineWidth', 2);
hold on;
plot(t, y2, 'b', 'LineWidth', 2);
plot(t, z2, 'g', 'LineWidth', 2);
xlabel('Time (Hours)');
ylabel('Concentration (arbitrary unit)');
legend('x (mRNA Concentration)','y (Protein Concentration)','z (Inhibiting Protein Concentration)');
title('Time course of x, y and z in presence of perturbation');
xlim([0 100]);  
ylim([0 0.7]);
perturbation_time = 50; 
xline(perturbation_time, 'k--', 'LineWidth', 2, 'DisplayName', 'Perturbation')
hold off;

%3D Phase Portrait
figure(2);
hold on;

plot3(x2(1:ceil(50 / dt)), y2(1:ceil(50 / dt)), z2(1:ceil(50 / dt)), 'k', 'LineWidth', 2);
plot3(x2(ceil(50 / dt):end), y2(ceil(50 / dt):end), z2(ceil(50 / dt):end), 'c', 'LineWidth', 2);
xlabel('x (mRNA concentration)');
ylabel('y (Protein concentration)');
zlabel('z (Inhibiting protein concentration)');
title('3D trajectory of system');
legend('Trajectory before perturbation', 'Trajectory after perturbation');
hold off;
grid on;

%% Simulating the Dynamics of the 2-element system
%% Define the parameters

%Definte time vector and steps
K = 0.1;    % Threshold of sigmoid
n = 10;    % Hill coefficient
dt = 1 / 3600;   % Time step in seconds
t = 0:dt:(100 - dt); %Time vector in hours
num_steps = length(t); % Number of steps

% Initialize vectors 
x = zeros(1, num_steps);
z = zeros(1, num_steps);

% Initial conditions 
x(1) = 0;
z(1) = 0;

%% Define the ODE system using Euler's method
%Integration without perturbation
for i = 1:num_steps-1
    dxdt = ((K^n) / (K^n + z(i)^n)) - x(i);
    dzdt = x(i) - z(i);
    
    x(i+1) = x(i) + dxdt * dt;
    z(i+1) = z(i) + dzdt * dt;
end


%% Plot the results

%Plot for 2D dynamics 
figure (3);
plot(t, x, 'r', 'LineWidth', 2);
hold on;
plot(t, z, 'g', 'LineWidth', 2);
xlabel('Time (Hours)');
ylabel('Concentration (arbitrary unit)');
legend('x (mRNA Concentration)', 'z (Inhibiting Protein Concentration)');
title('Time course of x and z in the 2-element system');
xlim([0 10]);  
ylim([0 0.7]);
hold off;

figure(4);  
plot(z, x, 'b', 'LineWidth', 2);  
xlabel('z (Inhibiting Protein Concentration)'); 
ylabel('x (mRNA Concentration)'); 
title('2D Phase Trajectory of x and z');  
grid on;  

%% Short Answer Questions

% 1. Equation choice B with K = 0.1 and n = 10 are the options that lead to
% stable oscillations in the system

% 2. As mentioned above, the above options result in stable oscillations
% which means after perturbation the system will return back to its
% previous oscillatory behavior. More specifically, at 50 hours the
% concentration of mRNA(x), the concentration of protein(y), and the
% concentration of inhibiting protein(z) are forced to a perturbation value
% of 0.13. Looking closely at the plot, the concentration of mRNA(x)
% increases the most and the concentration of protein(y) increases
% slightly. On the other hand, the concentration of inhibiting protein(z)
% decreases. As the system recovers you can see as the concentration of
% mRNA(x) increases so does the concentration of protein(y) and eventually
% the concentration of inhibiting protein(z) that results in the
% concentration of mRNA(x) decreasing. This is a textbook negative feedback
% loop with the purpose of bringing the system back to normal. As can
% be seen in Figure 1(plot 2) the negative feedback loop results in the
% system returning back to its normal oscillatory motion.

% 3. The dynamics between the 2-element system and the 3-element system are
% vastly different. The main reason for this is that there is no
% intermediate protein(protein y) in the 2 element system. Without the
% intermediate protein, the mRNA(x) is directly translating the inhibitory
% protein(z). This means that the gene responsible for transcribing
% mRNA(x) and translating inhibitory protein(z) has direct negative feedback
% on itself. Therefore, as the concentration of mRNA(x) increases the
% concentration of inhibitory protein(z) directly increases, subsequently,
% decreasing the concentration of mRNA(x) which directly decreases the
% concentration of inhibitory protein(z). Without an intermediate
% protein(y), there is no buffer that allows the mRNA concentration(x) to
% increase for a period of time without directly causing the inhibition of
% itself. As a result, the 2-element system has unstable dynamics as can be
% seen in Figure 3 where the concentration of mRNA(x) and inhibitory
% protein(z) are approaching very close to equilibrium around a
% concentration of ~0.12. To clarify, I use the word equilibrium as the
% concentrations of x and z are close to constant(0.12) and dx/dt and dz/dt
% are equal to 0 as 't' goes to infinity. In conclusion the unstable
% oscillations, or lack of oscillations, in the 2-element model vastly
% differ from the stable oscillations seen in the 3-element model.