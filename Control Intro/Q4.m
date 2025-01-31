%% Question 4

%% Question 4a

% Define the parameters of the system
run('parameters.m');

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

