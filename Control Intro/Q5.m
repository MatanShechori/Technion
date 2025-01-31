%% Question 5

run('parameters.m');
% Define the PI controller
Ki = [1,-1,10,0]; % Integral gain (you can adjust this value)
for i = 1:length(Ki)
    PI_controller = 1 + Ki(i)/s;
    P_PI = PI_controller * P;
    PM_PI = PI_controller * PM;

    %calculate poles and zeros for P_PI
    pp = pole(P_PI);
    zz = zero(P_PI);
    disp(['poles are at:    ', num2str(pp')]);
    disp(['zeros are at:    ', num2str(zz')]);

    figure;
    subplot(1,2,1);
    rlocus(P_PI);
    title(['PI Controller positive $k_{p}$ , $k_{i}$ = ', num2str(Ki(i))],'Interpreter','latex');

    subplot(1,2,2);
    rlocus(PM_PI);
    title(['PI Controller negative $k_{p}$ , $k_{i}$ = ', num2str(Ki(i))],'Interpreter','latex');
end

%% Question 5d

kp = 0.35;
ki = 0.01;
PI_controller = kp*(1 + ki/s);
P_PI=PI_controller*P;
T = P_PI/(1+P_PI);
figure;
step(T);