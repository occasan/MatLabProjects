% MCDB 108C, Spring 2024
% Coding assignment #5: Simulating a variant of Luria-Delbruck in-silico experiment
% Anthony Sacco

% alpha represents the probability of acquiring resistance at each cell
% Run the simulations with value: alpha = 1e-9; 
alpha = 1e-9;
beta = 1;

% parameter values
number_trials = 10000; % Number of experiments analyzed. Each experiment consists in the virtual counting of the number of resistant colonies.
number_generations = 21; % Number of generations of bacteria before the challenge is applied and the number of resistant colonies is counted. 

% initial values
number_resistant_colonies_Darwin = zeros(number_trials,1); % Number of resistant colonies observed for Darwinian model.
number_resistant_colonies_Lamarck = zeros(number_trials,1); % Number of resistant colonies observed for Lamarckian model.

% Compute the total number of bacteria expected after number_generations
% for Lamarckian model
final_number_bacteria = 150*2^21;

for j = 1:number_trials % Loop on the total number of experiments.
 
NWild = zeros(1,number_generations+1);% Vector containing the number of wild-type (nonresistant) bacteria at each generation of a given experiment.
NMutant = zeros(1,number_generations+1);% Vector containing the number of resistant type bacteria at each generation.
NWild(1) = 150; %Initial value number of wild-type
NMutant(1) = 0; %Initial value number of mutants

for i=1:number_generations % Loop on the number of generation of a given experiment.
    New_born_mutants = poissrnd(alpha*NWild(i)); % Use the Poisson disribution to simulate the number of resitant mutants that appeared during one generation change.
    New_lost_mutants= binornd(New_born_mutants, beta);
    NWild(i+1)=(2*NWild(i)) - New_born_mutants + New_lost_mutants; % Replication and removal of new born resistant bacteria.
    NMutant(i+1)=(2* NMutant(i)) + New_born_mutants - New_lost_mutants; % Update the number of resistant bacteria. 
end

% Update the number of resistant colonies observed for a given
% experiment for the Darwinian model.
number_resistant_colonies_Darwin(j)= NMutant(end);
% Update the number of resistant colonies observed for a given
% experiment for the Lamarckian model.
number_resistant_colonies_Lamarck(j) = poissrnd(alpha*final_number_bacteria); % For the Lamarckian model, the assumption is that resistance only appears during the challenge at the end of the experiment. 

end

% Plot and compare the results of both models.
close all
figure
subplot(2,1,1)
histogram(number_resistant_colonies_Darwin,-0.5:1:25.5,'Normalization','probability') % Plot the histrogram with the resuls of the Darwinian model
xlabel('Number of resistant colonies')
ylabel('Percentage (probability)')
title('Results of Darwinian model')
subplot(2,1,2)
histogram(number_resistant_colonies_Lamarck,-0.5:1:25.5,'Normalization','probability') % plot the histrogram with the resuls of the Lamarckian model
xlabel('Number of resistant colonies')
ylabel('Percentage  (probability)')
title('Results of Lamarckian model')

%% --> (HW5) Write your answers to the following questions below
% What do you notice in your histograms upon changing beta to equal, 1e-3, .1, and 1? 
% For each new adjustment of your beta variable, 
% how would you explain the changes or lack thereof that you notice between the results you get? 
% Why were the poisson and binomial distributions used for new_born_mutants and new_lost_mutants respectively?

% There were only small changes when increasing the beta-values
% respectively. Furthermore, these changes where only observed in the
% "jackpot" portion of the histogram and also were only observed on the
% Darwinian model. The changes noticed were a decreasing amount of bars on
% the right hand portion of the graph, meaning there was less probability
% for more resistant colonies. The most notable of these changes was when
% beta was increased from 0.1 to 1. These changes can be explained by both
% the beta value and what is happening in the for loop. When the
% probability of losing a drug resistant mutation(beta) is increased to a
% value such as 1 an increasing amount of New_born_mutants are lost to
% New_lost_mutants thus the mutation is not able to reproduce. This can
% clearly be seen in the histogram as a higher number of drug resistant
% colonies are decreasing with an increasing beta value. Similarily, the
% opposite can be explained when there is a lower beta value. This means
% less New_born_mutants are being lost and thus are able to reproduce which
% is present in the histogram as there is a greater number of drug
% resistant colonies

% A poisson distribution was used for New_born_mutants because there was a
% very large number of trials meaning a lot of bacteria were able to
% reproduce and there was a very low probability of acquiring a
% drug-resistant mutation. Both of these conditions being met allows for
% the use of a poisson distribution that is more accurate. A binomial
% distribution was used for New_lost_mutants because there is a low number
% of trials and a higher probability of losing a drug-resistant mutation
% with an increasing beta value. Since the condition of a high number of
% trials with a low amount of "events" for a poisson distribution is not
% met, it makes more sense to use a binomial distribution for more accurate
% results. 