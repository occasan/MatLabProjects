function [time_log,molecule_counter]=BirthDeathProcess_continuous(max_time,prod_rate,deg_rate)
% Script written as part of 108C class
% This function simulates  a birth and death process corresponding to the reaction: 
% gene -> mRNA -> deradation 
% The synthesis and degradation are assumed to take place on a continous time scale with the wait times distributed as negative exponentials.   
% This function is intended to be called in a loop by the wrapper BirthDeathProcess_averaging
% Date: Spring 2018
% Matthieu Louis, mlouis@ucsb.edu
% With the assistance of Tie Bo Wu

if nargin == 0 % Default values. Test whether max_time is defined externally. nargin refers to the number of argument of the function (i.e., the number of parameters entered as an input)
    max_time=1e3;
    prod_rate=1; % Production rate (molecules/min)
    deg_rate=0.02; % Degradation rate (1/min)
end
    
% Initialization of parameter
current_time=0;
step_counter=1;
time_log=0;
molecule_counter=0;

while current_time < max_time % Use a while loop to increase the current time until it exceeds the maximum value max-time.
    
    wait_time=-log(rand)/(prod_rate+(deg_rate*molecule_counter(step_counter))); % Simulated wait time: makes use of the fact that production and degradation are two independent Poisson processes.
    prob_prod=prod_rate/(prod_rate+deg_rate*molecule_counter(step_counter)); % Propensity of mRNA production
    
    current_time=current_time+wait_time; % Update current time
    step_counter=step_counter+1; % Update the number of steps (reactions) associated with the experiment.
    time_log(step_counter)=current_time; % Add the current time to the time log

    if rand < prob_prod % Defines whether production takes place based on Monte Carlo method.
        molecule_counter(step_counter)=molecule_counter(step_counter-1)+1; % Implements production
    else
        molecule_counter(step_counter)=molecule_counter(step_counter-1)-1; % Implements degradation
    end
      
    
end
  
if nargin == 0 % Plot results only when the number of input argument into the function is 0 (i.e., when the function is executed directly from the editor)
    close all
    subplot(2,1,1)
    plot(time_log(1:500),molecule_counter(1:500)) % Plots only the first 500 transitions
    ylabel('Number of mRNA molecules')
    xlabel('Time (min)')
    title('Full time course of the number of mRNA molecules')
    subplot(2,1,2)
    plot(time_log(1:50),molecule_counter(1:50),'*-') % Plots only the first 50 transitions
    hold on
    plot(time_log(1:50),ones(50,1),'r*')
    legend('Number of mRNA molecules','Timestamp of reaction event','Location','northwest')
    title('Initial change in the number of mRNA molecules')
    ylabel('Number of mRNA molecules')
    xlabel('Time (min)')
    hold off
end