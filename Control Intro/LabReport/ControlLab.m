%% Intro to Conrol lap report
% This report is based on the project and parameters.m


matlist12 = {'q12 0.1', 'q12 0.2', 'q12 0.3', 'q12 0.05', 'q12 0.36'};
simmatlist12 = {'q12 0.1 sim', 'q12 0.2 sim', 'q12 0.3 sim', 'q12 005 sim', 'q12 0.36 sim'};
simmatlist18 = {'q18 0.1 sim', 'q18 0.4 sim', 'q18 0.8 sim', 'q18 1.1 sim'};
matlist18 = {'q18 0.1','q18 0.4','q18 0.8',  'q18 1.1'};
for i= 1:length(matlist12)
    
    load(matlist12{i});
    load(simmatlist12{i});
    figure;
    plot(out.y_sim.Time,out.y_sim.Data);
    hold on;
    plot(y.time,y.signals.values);
    yline(30, '--r', 'Reference'); % Add horizontal line at y=30
    hold off;
    title(matlist12{i});
    xlabel('Time [s]');
    ylabel('Amplitude');
    legend('Simulation','Experiment', 'Reference');
    grid on;
end
% plotting the step response of the system for PI controller
for i= 1:length(matlist18)
    load(matlist18{i});
    figure;
    plot(out.y_sim.Time,out.y_sim.Data);
    hold on;
    plot(y.time,y.signals.values);
    yline(30, '--r', 'Reference'); % Add horizontal line at y=30
    hold off;
    title(matlist18{i});
    xlabel('Time [s]');
    ylabel('Amplitude');
    legend('Input', 'Output', 'Reference');
end

sisotool(kp*P);