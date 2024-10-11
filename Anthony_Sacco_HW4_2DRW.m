% MCDB 108C, Spring 2024
% Coding assignment #4: Simulating 2-D Random Walks
% INSERT: Anthony Sacco

%% Initialize random walk parameters
number_steps = 400; % Number of individual displacements of each particle (1 step/min)
number_particles = 1000; % Number of virus released at time 0 from position (0,0)

pos_x = zeros(number_steps,number_particles); % INSERT: Stores x-position of all the particles for each time step
pos_y = zeros(number_steps,number_particles); % INSERT: Stores y-position of all the particles for each time step
pos_dto = zeros(number_steps,number_particles); % Empty variable used to store each particle's distance to the origin (dto) 

distance_of_barrier = 50; % variable defining the distance to the barrier
number_particles_outside = 0; % variable to record the number of particles that have crossed the barrier

% tendency to move right
drift = 0.5; % you will vary this line to see the effect of vented airflow to the right

% probabilities of directional movement
prob_right = 0.5 + drift; %INSERT: The probability is equal to 0.5 without the drift because each particle has an equal chance of moving left or right
prob_up = 0.5;

%% Execute random walk trajectories for all particles, and compute distance to origin
for k = 1:number_particles % Outer for loop
    for i = 2:number_steps % Inner for loop 

        %Random choice for the x-direction
        if rand < prob_right
            pos_x(i, k) = pos_x(i-1, k) + 1;  
        else
            pos_x(i, k) = pos_x(i-1, k) - 1;
        end

        %Random Choice for y-direction
        if rand < prob_up
            pos_y(i, k) = pos_y(i-1, k) + 1;  
        else
            pos_y(i, k) = pos_y(i-1, k) - 1;  
        end
    end 
end
distances = sqrt(pos_x(:, end).^2 + pos_y(:, end).^2);
%% Compute number of particles past the barrier at the final time step
for k = 1:number_particles
    outside_flag = 0; % Binary variable used to record whether a particle crosses the barrier
    for i = 1:number_steps

        if abs(pos_x(i, k)) > distance_of_barrier
            outside_flag = 1;
            break; 
        end
    end
    if outside_flag == 1
        number_particles_outside = number_particles_outside + 1;
    end
end
%% Plot initial, midpoint, and final positions and display proportion of particles that leave barrier

% Plot positions at t = 0, t = 200, t = 400
figure
hold on
scatter(pos_x(1,:), pos_y(1,:),'o','filled','r') % Plot initial positions of all particles in red
scatter(pos_x(200,:), pos_y(200,:),'o','filled','g') % Plot positions of all particles at the midpoint (t = 200) in green
scatter(pos_x(end,:), pos_y(end,:),'o','filled','b') % Plot final positions of all particles at t = 400 in blue
legend('Initial Position', 'Position at t = 200', 'Final Position')
title('Particle Trajectory For Time Steps')
xlabel('X-Position')
ylabel('Y-Position')

% Display the number of particles that leave the barrier
disp(['Proportion of particles outside the barrier: ', num2str(number_particles_outside/number_particles)])

%% For an additional challenge, uncomment the code below, and plot probability distributions of final positions in x and y dimensions
% figure
% subplot(1,2,1)
% histogram(%INSERT: correct variable and indexing,50,'Normalization','pdf') % INSERT: comment for this line of code
% title('Probability Distribution in x-dimension')
% xlabel(% INSERT: Correct axis title)
% ylabel(%INSERT: Correct axis title)
% subplot(1,2,2)
% histogram(%INSERT: correct variable and indexing,50,'Normalization','pdf') % INSERT: comment for this line of code
% title('Probability Distribution in y-dimension')
% xlabel(%INSERT: Correct axis title)
% ylabel(%INSERT: Correct axis title)

%%%

% INSERT: Using comments, explain the outcome of your simulation results. How much did a change in airflow 
% (increase in drift) change the midpoint and final position of viral
% particles compared to no airflow (drift = 0)? What difference did you
% notice in how many virus particles crossed the barrier?

% When there is an increase in drift(0, 0.05, 0.5) there is a major change in the midpoint and final 
% positions. More specifically, the increase in airflow resulted in a
% positive shift towards the common area and ultimately the barrier. The
% difference in the amount of virus particles that crossed the barrier
% increased by ~0.015, ~0.35, 1 respective to the drift amount. This
% suggests that increasing drift(air flow) caused an increase in the amount
% of virus particles that entered the uninfected area. More interestingly,
% an increase of just drift = 0.05(5%) caused an increase of virus
% particles that crossed the barrier by approximately 35%