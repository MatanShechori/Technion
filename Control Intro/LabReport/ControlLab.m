%% Intro to Conrol lap report
% This report is based on the project and parameters.m

run parameters.m
kp = 6;

% Simulate step response
t = 0:0.01:15; % Time vector from 0 to 5 seconds
[y, t] = step(kp*P, t);
figure;
plot(t,y);

% Plot the step response
figure;
plot(t, y, 'LineWidth', 1.5);
grid on;
title('Step Response with Kp = 0.0485');
xlabel('Time [s])');
ylabel('Amplitude');

matlist12 = {'q12 0.1', 'q12 0.2', 'q12 0.3', 'q12 0.05', 'q12 0.36'};
matlist18 = {'q18 0.1','q18 0.4','q18 0.8',  'q18 1.1'};
for i= 1:length(matlist12)
    
    load(matlist12{i});
    figure;
    
    plot(u.time,u.signals.values);
    hold on;
    plot(y.time,y.signals.values);
    hold off;
    title(matlist12{i});
    xlabel('Time [s]');
    ylabel('Amplitude');
end
for i= 1:length(matlist18)
    load(matlist18{i});
    figure;
    plot(u.time,u.signals.values);
    hold on;
    plot(y.time,y.signals.values);
    hold off;
    title(matlist18{i});
    xlabel('Time [s]');
    ylabel('Amplitude');
end

sisotool(kp*P);