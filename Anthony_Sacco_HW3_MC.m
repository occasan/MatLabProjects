% MCDB 108C, Spring 2024
% Coding assignment #3: Monte-Carlo Simulations
% INSERT: Anthony Sacco

%% Parameter set 
p_infection = 1e-4; % probability that a cell gets infected by a given virus
num_viruses = 1e3; % total number of viruses. 
num_trials = 1e4; % total number of cells 

%% Create Monte Carlo simulation algorithm 
infection_tracker = zeros(num_trials,1); % vector recording the number of 
% infections for each cell
for i = 1:num_trials
    for j = 1:num_viruses
        if rand < p_infection
            infection_tracker(i) = infection_tracker(i)+ 1;
        end
    end
end

%% Analyze results
infection_stats = zeros(4,1); % vector recording the number of cells with 0, 1, 2 and >2 infections
for i = 1:num_trials 
    if infection_tracker(i) == 0
        infection_stats(1) = infection_stats(1) + 1;
    elseif infection_tracker(i) == 1
        infection_stats(2) = infection_stats(2) + 1;
    elseif infection_tracker(i) == 2
        infection_stats(3) = infection_stats(3) + 1;
    elseif infection_tracker(i) > 2
        infection_stats(4) = infection_stats(4) + 1;
    end
end 

%% Plot results of infected cells
bar(0:3, infection_stats / num_trials, 'FaceColor', [0.2 0.2 0.5]);
xticks(0:3);
xticklabels({'0 infections', '1 infection', '2 infections', '>2 infections'});
title('Distribution of Viral Infections per Cell');
xlabel('Number of Infections');
ylabel('Proportion of Cells');
grid on;


%INSERT: Does your result make sense? Explain in a few sentences what you
%see and why. Comment on the shape of the distribution you obtain. Does it 
%look like a bell-shaped (Gaussian) curve? Another type of distribution?

%This bar graph shows the proportion of cells that are associated with the
%states we defined(0 infections, 1 infection, etc..). This result of the
%bar graph makes sense given our original parameters. More specifically, we
%have a low probability of infection compared to a much higher number of
%cells which serves as a model to show how unlikely multiple infection
%events are. The bar graph would associate more with a Poisson distribution
%rather than a Gaussian distribution as most cells are observed with 0
%infections and the number of cells decreases as the infection count
%increases.

%% Excerise 1 Week 4 Discussion
clearvars

num_gfps = 3;
counter = zeros(num_gfps + 1, 1);

for i = 0: num_gfps
    counter(i+1) = nchoosek(num_gfps, i);
end

bar(counter);
xlabel = ('Number of Gfps in Daughter Cell A');
ylabel = ('Occurences');

%% Excersise 2 Week 4 Discussion
clearvars

% Parameters
n = 10; % Total GFPs
k = 3;  % GFPs to Daughter A
p = 0.5; % Probability of each GFP going to Daughter A

% Binomial distribution formula
P_k = nchoosek(n, k) * p^k * (1-p)^(n-k);

% Display the result
fprintf('The probability of Daughter A inheriting exactly 3 GFPs is %.4f\n', P_k);

