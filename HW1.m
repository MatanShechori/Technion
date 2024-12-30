%% Known Var of the Question

m = 1; % kg
g = 9.81; %m/sec^2
a = 1; % m
k = 4; %N/m
c1 = 0.1;
c3 = 0.001;
teta = 30*pi/180; % rad

%% quesion 5
L0 = 0.2;
f = @(x,L) k*(1- L/(sqrt(x^2+a^2)))*x-m*g*sin(teta);
initial_guess = -10:0.1:10;
root = zeros(size(initial_guess));
for i = 1:length(initial_guess)
    root(i)  = fzero(@(x)f(x,L0),initial_guess(i));
end
root = unique(round(root,5));

% Check stability of the solutions: 
q5_stable = is_stable(f,root,L0);
disp(q5_stable);


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

% find diffferent solution quantity based on L0
x = -10:0.01:10;
% isolating L0 from f when f=0
L0 = @(x) sqrt(a^2 + x^2) * (1 - (m*g*sin(teta))/(k*x));
L0_values = double(zeros(size(x)));

%Calculate L0 for each x
for i = 1:length(x)
    L0_values(i) = L0(x(i));
end

%plot L0 vs x
figure
plot(x,L0_values)
xlabel('x')
ylabel('L0')
title('L0 vs x')
ylim([0 20])
xlim([-10 10])
grid on
grid minor

% second plot for V(x) 
L0_cr = 3.1429;
x = -9:0.01:11;
L0 = [0.2,2.5,L0_cr,5,7];
V = @(x,L) 0.5*k*(sqrt(a^2 + x.^2) - L).^2 - m*g*sin(teta)*x;
tempV = zeros(size(x));
figure
for i = 1:length(L0)
    for j = 1:length(x)
        tempV(j) = V(x(j),L0(i));
    end
    plot(x,tempV)
    hold on
end

xlabel('x[m]','Interpreter','latex','FontSize',14);
ylabel('V(x)[j]','Interpreter','latex','FontSize',14);
legend('$\ell_{0}$ = 0.2','$\ell_{0}$ = 2.5','$\ell_{0}$ = 3.1429','$\ell_{0}$ = 5','$\ell_{0}$ = 7','Interpreter','latex','FontSize',8);
grid on
grid minor







%% Question 8


%% Question 9


%% Question 10

%% Question 11


%% Question 12


%% Functions

function [stable] = is_stable(f,test_points,L0)
    epsilon = 1e-6; % Small perturbation for numerical derivative
    stable =[test_points;zeros(size(test_points))];
    for i = 1:length(test_points)
    x = test_points(i);
    f_prime = (f(x + epsilon, L0) - f(x - epsilon, L0)) / (2 * epsilon); % Numerical derivative
        if f_prime < 0
        stable(2,i) = 0;
        else
        stable(2,i) = 1;
        end
    end
end