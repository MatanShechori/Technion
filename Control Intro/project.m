%% Question 4

%% Question 4a

% Define the parameters of the system
run('parameters.m');

% Create a figure with 2 subplots for P and -P
figure;

% First subplot for P
subplot(1, 2, 1);
rlocus(P);
title('Root Locus for positive Kp');

% Second subplot for -P
figure;
subplot(1, 2, 2);
rlocus(-P);
title('Root Locus for negative Kp');

pol = pole(P);
zer = zero(P);

disp('Poles of the system are:');
disp(pol);
disp('Zeros of the system are:');
disp(zer);


%% Question 5

% set the controller to be a PI controller

Kp =0.0485;
Ki = [2,4];
for i = 1:length(Ki)
    figure;
    C = Kp*(1+(Ki(i)/s));
    figure;
    subplot(1,2,1)
    rlocus(Cp*P);
    title(['Root Locus for Ki = ', num2str(Ki(i))]);
    subplot(1,2,2)
    rlocus(-Cp*P);
    title('root locus for negative Kp');
end

