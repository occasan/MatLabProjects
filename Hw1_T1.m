% MCDB 108C, Spring 2024
% Coding assignment #1
% Anthony Sacco - Solution for beginner track

clearvars; 

dx = 0.25; 
x_values = -5 : dx : (5 - dx);

mu = 0;
sigma = 1;

sum_total = 0.;
for i = 1:length(x_values)
    x = x_values(i);
    partial_sum = (1 / sqrt(2 * pi)) * exp(-x^2 / 2);
    sum_total = sum_total + partial_sum * dx; 
end

disp('---------------------------------------------------')
disp('The value of the sum is:')
disp(sum_total)