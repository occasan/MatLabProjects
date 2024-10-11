% If, else, and else if statements in matlab
randnum = randi(4);
if randnum == 1
    disp('The number is 1')
elseif randum ^= 1
    disp('Something else not 1')
else
end

%% Excerise 1
clearvars;

for i = 1:20
r = randi([1, 10]);

    if r <= 7
        disp('No luck today')
    elseif r>7 r<10;
        disp('Good luck today')
    else 
        disp('Extremely goodluck today')
    end
end
%% Excerise 2

a = [111, 33, -1, 2 ,-5, 3, 4, 10, -22];

p_array = [];
n_array = [];

p_i = 1;
n_i = 1;

for i = 1:length(a)
    if a(i) < 0
        n_array(n_i) = a(i);
        n_i = n_i+1;
    else
        p_array(p_i) = a(i)
        p_i = p_i+1;
    end
end

disp('Negative values')
disp(n_array)

disp('Positive values')
disp(p_array)
        


