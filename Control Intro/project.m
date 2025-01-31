%% Question 4

%% Question 4a

% Define the parameters of the system
run('parameters.m');

s = tf('s');
% Define the transfer function of the system
PM = -P;
pp = pole(P);
disp(['poles are at:    ', num2str(pp')]);

% Create root locus plots for P and PM
figure;
subplot(1,2,1);
rlocus(P);
title('$k_{p} > 0$','Interpreter','latex');

subplot(1,2,2);
rlocus(PM);
title('$k_{p} < 0$','Interpreter','latex');
%% Question 5

run('parameters.m');
% Define the PI controller
Ki = [1,-1,10,0]; % Integral gain (you can adjust this value)
for i = 1:length(Ki)
    PI_controller = 1 + Ki(i)/s;
    P_PI = PI_controller * P;
    PM_PI = PI_controller * PM;

    figure;
    subplot(1,2,1);
    rlocus(P_PI);
    title(['PI Controller positive $k_{p}$ , $k_{i}$ = ', num2str(Ki(i))],'Interpreter','latex');

    subplot(1,2,2);
    rlocus(PM_PI);
    title(['PI Controller negative $k_{p}$ , $k_{i}$ = ', num2str(Ki(i))],'Interpreter','latex');
end
