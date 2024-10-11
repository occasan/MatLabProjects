%% Coding Assignment #7, MCDB 108C, Spring 2024 
% Description: 1-ELEMENT POSITIVE FEEDBACK LOOP
% Anthony Sacco
%% Parameters
max_time = 20; % In minute
dt = 0.01;
time_vect = 0:dt:max_time;
n_step = length(time_vect);

gamma = 20; % Production rate (molecules/min) 
theta = 10; % Threshold of sigmoid
k = 1; % Degradation rate (1/minute)
n_hill = 10; % Hill coefficient of sigmoid

initial_x = 12; % Initial condition for deterministic and stochastic model

%% Part 1: Deterministic ODE model
x = zeros(n_step, 1); % Initialize output of deterministic model (Euler's integration)
x(1) = initial_x; % Initial condition x(0)

for i = 1:(n_step-1)

    %dxdt = gamma - k*x(i);
    dxdt = gamma * (x(i)^n_hill)/((theta^n_hill)+x(i)^n_hill - k*x(i));
    x(i+1) = x(i) + dxdt*dt;    

end

%% PART 2: Stochastic birth and death model
% The code below is based on Script 4: Continuous BDP seen in the class
current_time = 0;
step_counter = 1;
time_log = 0;
molecule_counter = initial_x;

while current_time < max_time % Use a while loop to increase the current time until it exceeds the maximum value max-time.
    
    current_prod_rate = gamma; %INSERT 2
    %current_prod_rate = gamma * (molecule_counter(step_counter)^n_hill)/((theta^n_hill)+molecule_counter(step_counter)^n_hill);
    current_deg_rate = k*molecule_counter(step_counter); %INSERT 3
    %current_deg_rate = k*molecule_counter(step_counter);

    wait_time=-log(rand)/(current_prod_rate+current_deg_rate); % Simulated wait time: makes use of the fact that production and degradation are two independent Poisson processes.
    prob_prod=current_prod_rate/(current_prod_rate+current_deg_rate); % Propensity of production reaction
    
    current_time=current_time+wait_time; % Update current time
    step_counter=step_counter+1; % Update the number of steps (reactions) associated with the experiment.
    time_log(step_counter)=current_time; % Add the current time to the time log

    if rand < prob_prod % Defines whether production takes place based on Monte Carlo method.
        molecule_counter(step_counter)=molecule_counter(step_counter-1)+1; % Implements production
    else
        molecule_counter(step_counter)=molecule_counter(step_counter-1)-1; % Implements degradation
    end   
    
end

%% Part 3: Plotting and analyzing the results
clf
plot(time_vect, x, 'LineWidth', 2)
hold on
plot(time_log,molecule_counter,'r','LineWidth', 2)
xlabel('Time (min)', 'FontSize', 14)
ylabel('Number of molecules', 'FontSize', 14)
legend('Result of deterministic ODE model','Result of stochastic model')
xlim([0 20])
ylim([0 30])

%% COMMENT 1: Describe the behavior of the non-regulated thermostat for the deterministic and stochastic models

% For the deterministic model, the non-regulated thermostat shows
% a continuous increase in the number of molecules over time until it stabilizes 
% at a steady state. This steady state is reached when the production and degradation 
% rates are equal. In contrast, the stochastic model displays 
% significant fluctuations around the mean behavior. The number of molecules 
% randomly jumps and drops due to the probabilistic nature of the molecular 
% interactions, resulting in a noisy trajectory compared to the 
% smooth curve of the deterministic model.


%% COMMENT 2: Describe the behavior of the regulated positive feedback loop (PFL) for the deterministic and stochastic models

% In the regulated positive feedback loop (PFL), the deterministic model 
% initially shows a smooth increase in the number of molecules before stabilizing 
% at a higher concentration due to the feedback mechanism. The feedback introduces 
% non-linearity, leading to a stable state at a higher level than the initial 
% condition. The stochastic model with PFL exhibits pronounced fluctuations compared to before. 
% The feedback amplifies the randomness, causing variability and 
% occasional rapid changes in the number of molecules, making the system less 
% predictable compared to the deterministic model.

% Comment 3: Before my deterministic model for both cases represented a
% sigmoidal shape. However, after playing with the code some more both of
% my deterministic models are stuck being completely linear despite
% reverting the changes. Quite confused and was not able to figure this out
% before the assignment was due. Hopefully, I do not lose points for this as
% I was still able to answer the questions. Thank you!








