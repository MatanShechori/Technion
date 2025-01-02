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
title('Q7 figure 1: L0 vs x')
ylim([0 20])
xlim([-10 10])
grid on
grid minor

% second plot for V(x) 
L0_cr = 3.1429;
x = -9:0.01:11;
L0_samples = [0.2,2.5,L0_cr,5,7];
V = @(x,L) 0.5*k*(sqrt(a^2 + x.^2) - L).^2 - m*g*sin(teta)*x;
tempV = zeros(size(x));
figure
for i = 1:length(L0_samples)
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
title('Q7 figure 2: V(x) vs x');


%% Question 8

% Plot x vs L0 with dashed line after L0 = L0_cr
figure
x_cr = -(m*g*sin(teta)*a^2/k)^(1/3);
x1 = -10:0.01:x_cr;
L1_values = zeros(size(x1));
x2 = x_cr:0.01:10;
L2_values = zeros(size(x2));
for i = 1:length(x1)
    L1_values(i) = L0(x1(i));
end
plot(L1_values,x1, 'b')
xlim([-15 15])
ylim([-10 10])
hold on
for i = 1:length(x2)
    L2_values(i) = L0(x2(i));
end
plot(L2_values,x2,'b--')
grid on
grid minor
title('Q8 figure 1: x vs L0')


%% Question 9
% state space representation of the system
time_vec  = 0:0.01:25;
l0 =0.2;
x0 = [-2,-1,0,1,2];
xdot = @(t,x,l0) [x(2); g*sin(teta) - (c1*x(2)+c3*x(2)^3)/m - k*(1-l0/sqrt(x(1)^2+a^2))*x(1)/m];

% Figure 1 l0 = 0.2:
figure;
for i = 1:length(x0)
    [~,X] = ode45(@(t,x) xdot(t,x,l0),time_vec,[x0(i) 0]);
  
    plot(X(:,1),X(:,2));
    hold on
end
plot(q5_stable(1,1), 0, 'r.','MarkerSize', 20);
legend('$x_{0} = -2$','$x_{0} = -1$','$x_{0} = 0$','$x_{0} = 1$','$x_{0} = 2$','stable eq','Interpreter','latex','FontSize',8);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',14);
title('Q9: $l_{0} = 0.2$','Interpreter','latex','FontSize',14);
grid on
grid minor
xlim([-4 6])
ylim([-8 8])

% figure 2 l0 = 5:
    l0 = 5;
    figure;
for i = 1:length(x0)
     [~,X] = ode45(@(t,x) xdot(t,x,l0),time_vec,[x0(i) 0]);
  
      plot(X(:,1),X(:,2));
       hold on
end
plot(q6_stable(1,1), 0, 'r.','MarkerSize', 20)
plot(q6_stable(1,2), 0, 'bo')
plot(q6_stable(1,3), 0, 'r.','MarkerSize', 20)
legend('$x_{0} = -2$','$x_{0} = -1$','$x_{0} = 0$','$x_{0} = 1$','$x_{0} = 2$', 'stable eq','unstable eq', 'Interpreter','latex','FontSize',8);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',14);
title('Q9: $l_{0} = 5$','Interpreter','latex','FontSize',14);
xlim([-6 12])
ylim([-12 12])
grid on
grid minor


%% Question 10

% Initial conditions
x_eq = q6_stable(1,3);
l0 = 5;
k_eq = k*(1-l0/sqrt(x_eq^2+a^2))+k*l0*x_eq^2/(x_eq^2+a^2)^1.5;
Omegan = sqrt(k_eq/m);
zeta = c1/(2*m*Omegan);
Omegad = Omegan*sqrt(1-zeta^2);


u0 = 2;
x0 = u0+x_eq;
v0 = 0;
Amp = sqrt(u0^2+(v0+zeta*Omegan*u0)/Omegad^2);
phi0 = atan((u0*Omegad)/(v0+zeta*Omegan*u0));

time_vec = 0:0.01:100;

% Response for u0 = 2:
u = zeros(size(time_vec));
for i = 1:length(time_vec)
    u(i) = Amp*exp(-zeta*Omegan*time_vec(i))*sin(Omegad*time_vec(i)+phi0);
end

[~,x] = ode45(@(t,x) xdot(t,x,l0),time_vec,[x0 0]);
figure
plot(time_vec,x(:,1)+x_eq);
hold on
plot(time_vec,u+x_eq,'--g')
grid on
grid minor



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