%% Known Var of the Question

m = 1; % kg
g = 9.81; %m/sec^2
a = 1; % m
k = 4; %N/m
c1 = 0.1;
c3 = 0.001;
teta = 30*pi/180; % rad

%% quesion 5

f = @(x) k*(1- L0/(sqrt(x^2+a^2)))*x-m*g*sin(teta);

a = 1; % Lower bound
b = 3; % Upper bound
tol = 1e-6; % Tolerance for convergence

while (b - a) / 2 > tol
    c = (a + b) / 2; % Midpoint
    if f(c) == 0 % Check if we found the root
        break;
    elseif f(a) * f(c) < 0
        b = c; % Root is in the left half
    else
        a = c; % Root is in the right half
    end
end

root = (a + b) / 2;
disp(root);


%% question 6
% find all the available solutions of f(x) when L0 = 5 m
L0 = 5;
for i = 1:length(initial_guess)
    root(i)  = fzero(@(x)f(x,L0),initial_guess(i));
end
root = unique(round(root,5));

% Check stability of the solutions: 
q6_stable = is_stable(f,root,L0);
disp(q6_stable);


%% question 7



%% Q8


%% Q9


%% Q10

%% Q11


%% Q12

function [stable] = is_stable(f,test_points,L0)
    epsilon = 1e-6; % Small perturbation for numerical derivative
    stable =[test_points;zeros(size(test_points))];
    for i = 1:length(test_points)
    x = test_points(i);
    f_prime = (f(x + epsilon, L0) - f(x - epsilon, L0)) / (2 * epsilon); % Numerical derivative
        if f_prime <= 0
        stable(2,i) = 0;
        else
        stable(2,i) = 1;
        end
    end
end