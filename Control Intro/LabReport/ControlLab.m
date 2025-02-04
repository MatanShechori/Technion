

% Pre-load all data for q12
matlist12 = {'q12 0.05','q12 0.1','q12 0.2','q12 0.3',  'q12 0.36'};
simmatlist12 = {'q12 005 sim','q12 0.1 sim', 'q12 0.2 sim', 'q12 0.3 sim',  'q12 0.36 sim'};
simmatlist18 = {'q18 0.1 sim', 'q18 0.4 sim', 'q18 0.8 sim', 'q18 1.1 sim'};
matlist18 = {'q18 0.1','q18 0.4','q18 0.8', 'q18 1.1'};
for i = 1:length(matlist12)
    % Load experimental data
    temp = load(matlist12{i});
    exp_data(i).time = temp.y.time;
    exp_data(i).values = temp.y.signals.values;
    
    % Load simulation data
    temp = load(simmatlist12{i});
    sim_data(i).time = temp.out.y_sim.Time;
    sim_data(i).values = temp.out.y_sim.Data;
end

% Plot q12 data
for i = 1:length(matlist12)
    figure;
    plot(sim_data(i).time, sim_data(i).values);
    hold on;
    plot(exp_data(i).time, exp_data(i).values);
    yline(30, '--r', 'Reference');
    hold off;
    title(matlist12{i});
    xlabel('Time [s]');
    ylabel('Amplitude');
    legend('Simulation', 'Experiment', 'Reference');
    grid on;
end

% Pre-load all data for q18
for j = 1:length(matlist18)
    % Load experimental data
    temp = load(matlist18{j});
    exp_data18(j).time = temp.y.time;
    exp_data18(j).values = temp.y.signals.values;
    
    % Load simulation data
    temp = load(simmatlist18{j});
    sim_data18(j).time =  temp.out.y_sim.time;
    sim_data18(j).values = temp.out.y_sim.Data;
end

% Plot q18 data
for j = 1:length(matlist18)
    figure;
    plot(sim_data18(j).time, sim_data18(j).values);
    hold on;
    plot(exp_data18(j).time, exp_data18(j).values);
    yline(30, '--r', 'Reference');
    hold off;
    title(matlist18{j});
    xlabel('Time [s]');
    ylabel('Amplitude');
    legend('Simulation', 'Experiment', 'Reference');
    grid on;
end

