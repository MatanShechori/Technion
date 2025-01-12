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

% Define resolution and x range

x = -10:1e-2:10;

% Initialize L0 values
L0_values = double(zeros(1, length(x)));

% Calculate L0 for each x
for i = 1:length(x)
    L0_values(i) = sqrt(a^2 + x(i)^2) * (1 - m * g * sin(teta) / (k * x(i)));
end

% Define critical values
x_cr = -(m * g * sin(teta) * a^2 / k)^(1/3); % Critical x value
L0_cr = sqrt(a^2 + x_cr^2) * (1 - (m * g * sin(teta)) / (k * x_cr)); % Critical L0 value

% Find indices for critical and zero points
ind_cr = round((x_cr + 10) / 0.01);
ind_0 = find(x == 0);

% Separate stable and unstable branches
stable1 = [x(1:ind_cr-1); L0_values(1:ind_cr-1)];
stable2 = [x(ind_0+1:end); L0_values(ind_0+1:end)];
unstable = [x(ind_cr+1:ind_0-1); L0_values(ind_cr+1:ind_0-1)];

% Plot the results
figure;
plot(unstable(2,:), unstable(1,:), 'r--', 'LineWidth', 2); % Unstable branch
hold on;
plot(stable1(2,:), stable1(1,:), 'b', 'LineWidth', 1.5); % Stable branch (1st)
plot(stable2(2,:), stable2(1,:), 'g', 'LineWidth', 1.5); % Stable branch (2nd)
plot(L0_cr, x_cr, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Critical point

% Add labels, legend, and grid
xlabel('$L_0$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$x_{eq}$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
legend('Unstable points', 'Stable points (1st)', 'Stable points (2nd)', 'Critical point', ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);
grid on;
grid minor;
xlim([-15 15]);
ylim([-10 10]);
title('Q8: Equilibrium Points vs $L_0$', 'Interpreter', 'latex', 'FontSize', 14);


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
plot(time_vec,x(:,1));
hold on
plot(time_vec,u+x_eq,'--r')
grid on
grid minor
xlabel('Time [s]')
ylabel('x [m]')
title('Q10: Response for $u_{0} = 2$','Interpreter','latex','FontSize',14)

% Response for u0 = 6.7:
u0 = 6.7;
x0 = u0+x_eq;
Amp = sqrt(u0^2+(v0+zeta*Omegan*u0)/Omegad^2);
phi0 = atan((u0*Omegad)/(v0+zeta*Omegan*u0));

u = zeros(size(time_vec));
for i = 1:length(time_vec)
    u(i) = Amp*exp(-zeta*Omegan*time_vec(i))*sin(Omegad*time_vec(i)+phi0);
end

[~,x] = ode45(@(t,x) xdot(t,x,l0),time_vec,[x0 0]);
figure
plot(time_vec,x(:,1));
hold on
plot(time_vec,u+x_eq,'--r')
grid on
grid minor
xlabel('Time [s]')
ylabel('x [m]')
title('Q10: Response for $u_{0} = 6.7$','Interpreter','latex','FontSize',14)

%% Question 11
c = 1e-3;
Omegan = sqrt(k/m);
time_vec = 0:c:100;

%ode q11:
xdot11 = @(t,x) [x(2); -(k/m)*x(1)-(c/m)*abs(x(1))*x(2)^3];

%intial conditions:
x0 = [12;0];


[t1,y] = ode45(@(t,x)xdot11(time_vec,x),time_vec,x0);

figure
plot (t1,y(:,1));
xlabel('Time [s]')
ylabel('Amplitute [m]')
title('Q11: $x_{0} = 12$','Interpreter','latex','FontSize',14)
ylim([-12 12])
grid on
grid minor

%intial conditions:
x0 = [50,0];


[t1,y] = ode45(@(t,x)xdot11(time_vec,x),time_vec,x0);

figure
plot (t1,y(:,1));
xlabel('Time [s]')
ylabel('Amplitute [m]')
title('Q11: $x_{0} = 50$','Interpreter','latex','FontSize',14)
ylim([-60 60])
grid on
grid minor


%% Question 12

%part b
syms R real;
x_eq = q6_stable(1,3);
l0 = 5;
k_eq = k*(1-l0/sqrt(x_eq^2+a^2))+k*l0*x_eq^2/(x_eq^2+a^2)^1.5;
omegan_new = sqrt(k_eq/m);
w = linspace(0.1*Omegan,5*Omegan,100);
Omega = w/omegan_new;
M0 = 1:2:11;
response = zeros(length(M0),length(w));
figure
hold on
for i = 1:length(M0)
    F = a*M0(i)/(m*omegan_new^2*(x_eq^2+a^2));
    for j = 1:length(w)
        c_eq =c1+3*(c3/4)*(R*M0(i))^2*w(j)^2;
        delta = c_eq/(m*omegan_new);
        ratio = ((1-w(j)^2)^2+(delta*w(j))^2)*R^2-(F/M0(i))^2==0;
        sol = vpasolve(ratio,R,1);
        response(i,j) = double(abs(sol(1)));
    end
    plot(w,response(i,:)*M0(i));
end
xlabel('$\omega_{n}[rad/sec]$','Interpreter','latex','FontSize',14);
ylabel('$A[m]$','Interpreter','latex','FontSize',14);
set(gca, 'XScale', 'log'); 
legend(arrayfun(@(x) ['$M_0 = ', num2str(x), '$'], M0, 'UniformOutput', false), 'Interpreter', 'latex', 'FontSize', 8);
title('Q12: $A$ vs $\omega_{n}$','Interpreter','latex','FontSize',14);
grid on
grid minor

% part c
% simulate f response with added c_eq and M0=5 for specific frequencies and different initial conditions

M0 = 5;
frequencies = [0.1 * Omegan, 5 * Omegan];
initial_conditions = [0.1, 10];
res_t = 1e-2;

figure;
for f_idx = 1:length(frequencies)
    w = frequencies(f_idx);
    time_span = 0:res_t:300 * 2 * pi / w; % 100 cycles
    for ic_idx = 1:length(initial_conditions)
        x0 = initial_conditions(ic_idx);
        xdot12 = @(t, x) [x(2); (a * M0 * cos(w * t) / (a^2 + x_eq^2) - c1 * x(2) - c3 * x(2)^3 - k * x(1) * (1 - l0 / sqrt(a^2 + x_eq^2))) / m];
        [t, y] = ode45(xdot12, time_span, [x0; 0]);
        
        subplot(length(frequencies), 1, f_idx);
        plot(t, y(:, 1),'b');
        

        hold on;
        title(['$\omega = ', num2str(w), '$'], 'Interpreter', 'latex', 'FontSize', 10);
        xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 10);
        ylabel('x [m]', 'Interpreter', 'latex', 'FontSize', 10);
        ylim([-0.3 0.3]); % Adjust y-scale
        grid on;
        grid minor;
    end
end
legend(arrayfun(@(x) ['$x_0 = ', num2str(x), '$'], initial_conditions, 'UniformOutput', false), 'Interpreter', 'latex', 'FontSize', 8);

%% Question 12 d

% Ensure y is defined
if ~exist('y', 'var')
    error('Variable y is not defined. Ensure y is calculated before running this section.');
end

% Parameters
M_0 = 1; % Amplitude of the external force
omega = 0.1 * Omegan; % Example frequency, adjust as needed
t_end = 100; % Total simulation time (example value)
dt = 0.01; % Time step

% Time vector
t = 0:dt:t_end;

% External force
M = M_0 * cos(omega * t);

% Initial conditions
x0 = [0; 0]; % Initial position and velocity

% ODE function
odefun = @(t, x) [x(2);
    (M_0 * a / (x(1)^2 + a^2) - c1 * x(2) - c3 * x(2)^3 - k * x(1) * (1 - l0 / sqrt(x(1)^2 + a^2))) / m];

% Solve ODE
[t, x] = ode45(odefun, t, x0);

% Cut the last 10 revolutions
num_revolutions = 10;
T = 2 * pi / omega; % Period of one revolution
t_cut = t_end - num_revolutions * T;
indices = t >= t_cut;
tvec = t(indices);
x_in = M(indices)';
x_out = x(indices, 1);

% Ensure x_in and x_out are column vectors
x_in = x_in(:);
x_out = x_out(:);

% Model matrix
Model = [cos(omega * tvec), sin(omega * tvec), ones(size(tvec))];

% Solve least squares problem
consts = Model \ [x_in, x_out];
A_in = consts(1, 1) + 1i * consts(2, 1);
A_out = consts(1, 2) + 1i * consts(2, 2);

% Display results
disp(['Amplitude of input signal: ', num2str(abs(A_in))]);
disp(['Phase of input signal: ', num2str(angle(A_in))]);
disp(['Amplitude of output signal: ', num2str(abs(A_out))]);
disp(['Phase of output signal: ', num2str(angle(A_out))]);

%% Question 12 e

% Frequency range
omega = linspace(0.1, 10, 1000); % Adjust the range and resolution as needed

% Preallocate array for |A|
A_magnitude = zeros(size(omega));

% Define the constants for the equation
coeff = 3/4 * 10^(-3);

for idx = 1:length(omega)
    % Current frequency
    w = omega(idx);
    
    % Define the denominator of the equation
    term1 = (-w^2 / 2 + 0.382)^2;
    term2 = (w / 2)^2 * (coeff * w^2 - 0.1)^2;
    denominator = term1 + term2;
    
    % Calculate |A|
    A_magnitude(idx) = (0.0128 * M_0)^2 / sqrt(denominator);
end

% Create a figure with two subplots
figure;

% Subplot 1: Frequency response of the linearized system
subplot(1, 2, 1);
plot(omega, A_magnitude, 'LineWidth', 2);
xlabel('\omega (rad/s)');
ylabel('|A|[dB]');
title('Frequency Response of the Linearized System');
grid on;

% Subplot 2: Nonlinear system response
subplot(1, 2, 2);
plot(tvec, x_out, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Nonlinear System Response');
grid on;

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
