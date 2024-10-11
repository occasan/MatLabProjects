%% Coding Assignment #6, MCDB 108C, Spring 2024 
% Description: EXPLORING THE MERGING PROPERTY OF POISSON PROCESSES

%% Parameters
clearvars;
max_time = 120*60; % Total time steps (seconds) over which the production is simulated

%Define rates
slow_rate = 1; % rate constant of the production of a molecule, given in 1/second
fast_rate = 4;
merged_rate = slow_rate + fast_rate;

%% Part 1: Simulation of slow Poisson process

slow_wait_time_log = []; %Create an empty array to store wait times later
slow_total_time_log = []; %Create an empty array to store exact time later
%
slow_total_time = 0; %Starting time
i = 0; %Index

while slow_total_time < max_time 
    
    i = i+1;

    % Generate a random wait time for the slow Poisson process
    stochastic_variable = rand();
    slow_wait_time = -log(stochastic_variable)/slow_rate;
    slow_wait_time_log(i) = slow_wait_time;

    % Grow the vectors with data values
    slow_total_time = slow_total_time + slow_wait_time; %Add the wait time to the total amount of time
    slow_total_time_log(i) = slow_total_time; 
    
end

%% Part 2: Simulation of fast Poisson process 

fast_wait_time_log = [];
fast_total_time_log = [];
%
fast_total_time = 0;
i = 0;

while fast_total_time < max_time 
    
    i = i+1;

    % Generate a random wait time for the fast Poisson process
    stochastic_variable = rand(); 
    fast_wait_time = -log(stochastic_variable)/fast_rate; 
    fast_wait_time_log(i) = fast_wait_time; 

    % Grow the vectors with data values
    fast_total_time = fast_total_time + fast_wait_time;
    fast_total_time_log(i) = fast_total_time;
    
end

%% Part 3: Simulation of combined processes

slow_wait_time_log_merged = [];
slow_total_time_log_merged = [];
fast_wait_time_log_merged = [];
fast_total_time_log_merged = [];
%
merged_total_time = 0;
i_merged = 0;
i_fast = 0;
i_slow = 0;

fast_wait_time_merged = 0;
slow_wait_time_merged = 0;

while merged_total_time < max_time
    
    i_merged = i_merged+1;

    % Generate a random wait time and store value for the combined 
    stochastic_variable = rand();
    merged_wait_time = -log(stochastic_variable)/merged_rate;
    merged_total_time = merged_total_time + merged_wait_time;

    if rand < fast_rate/merged_rate
        i_fast= i_fast + 1;
        fast_wait_time_log_merged(i_fast) = fast_wait_time_merged + merged_wait_time;
        fast_total_time_log_merged(i_fast) = merged_total_time;
        fast_wait_time_merged = 0; %To ensure the next protein expression event is starting from time zero after the previous one
        slow_wait_time_merged = slow_wait_time_merged + merged_wait_time;
    else
        i_slow= i_slow + 1;
        slow_wait_time_log_merged(i_slow) = slow_wait_time_merged + merged_wait_time;
        slow_total_time_log_merged(i_slow) = merged_total_time;
        slow_wait_time_merged = 0; %To ensure the next protein expression event is starting from time zero after the previous one
        fast_wait_time_merged = fast_wait_time_merged + merged_wait_time;
    end
    
end

%% Part 4: Plotting and analyzing the results

%FIGURE 1
figure (1);
subplot(2,2,1);
histogram(fast_wait_time_log,0:0.2:5, 'Normalization', 'probability');
title('Fast production - individual');
ylabel ('Probability'); xlabel ('Wait times (second)');

subplot(2,2,3);
histogram(fast_wait_time_log_merged,0:0.2:5, 'Normalization', 'probability');
title('Fast production - combined');
ylabel ('Probability'); xlabel ('Wait times (second)');

subplot(2,2,2); %Plot the histograms for the individual slow rate vs. the merged slow rate
histogram(slow_wait_time_log,0:0.2:5, 'Normalization', 'probability') 
title('Slow production - individual');
ylabel ('Probability'); xlabel ('Wait times (second)');

subplot(2,2,4);
histogram(slow_wait_time_log_merged, 0:0.2:5, 'Normalization', 'probability'); 
title('Slow production - combined');
ylabel ('Probability'); xlabel ('Wait times (second)');

%FIGURE 2
% Plot the results of the individual simulations
figure (2);
max_index_fast_total_time=sum(fast_total_time_log<10);
scatter(fast_total_time_log(1:max_index_fast_total_time),1,'filled','or');

hold on;
max_index_slow_total_time=sum(slow_total_time_log<10);
scatter(slow_total_time_log(1:max_index_slow_total_time),1,'filled','ob');

% Plot the results of the merged simulations
max_index_fast_total_time_combined= sum(fast_total_time_log_merged < 10);
scatter(fast_total_time_log_merged(1:max_index_fast_total_time_combined),2,'*r');

max_index_slow_total_time_combined= sum(slow_total_time_log_merged < 10);
scatter(slow_total_time_log_merged(1:max_index_slow_total_time_combined),2,'*b');

ylim([0.8 2.2]);
xlabel('Time (s)');
title('o - single processes; x - combined processes','red - fast; blue - slow');

%% COMMENT: What is one advantage of simulating the merged Poisson process compared to simulating and then summing the individual Poisson processes?

% The main advtange of simulating the merged Poisson process is
% that by taking the sum of the fast rate and the slow rate to create the
% merged rate, you are effectively decreasing the amount of time between
% events. This means you see a greater amount of proteins being
% synthesized, whether fast or slow, in a given time period that can then
% be compared against eachother on the same visualization (ex. scatter
% plot). Ultimately, this is more efficient than trying to simulate two
% seperate Poisson processes and summing them together. 

