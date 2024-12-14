%% Known Var of the Question

m = 1; % kg
g = 9.81; %m/sec^2
L0 = 0.2; %m
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



%% question 7

%% Q8


%% Q9


%% Q10

%% Q11


%% Q12

