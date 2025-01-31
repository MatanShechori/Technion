%% Section 1: Thermal conduction:

% Create matlab tables:
plates = readtable('plates.csv');
fins = readtable('fins&ribs.csv');
plates = plates(~isnan(plates.Var1), :);
fins = fins(~isnan(fins.Var1), :);

% Question 1: Graph of temp vs distance from hotplate
average_temp = mean([plates.Li1, plates.Li2, plates.Li3, plates.Li4, plates.Li5], 2); % average temp from measurement
frame2dis = plates.Var1(end) / 207; % distance between frames [mm]
distance_vec = plates.Var1 ./ frame2dis;

figure
plot(distance_vec, average_temp, 'r')
xlabel('Distance from hotplate [mm]')
ylabel('Temperature [C]')
title('Temperature vs Distance from Hotplate')
grid on

steel = 8 * frame2dis;
granite = steel + 55;
aluminum = granite + 45;
glass = aluminum + 88;

% Create vertical lines on the figure to show the material change
hold on
xline(steel, '--b', 'Steel', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');
xline(granite, '--g', 'Granite', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');
xline(aluminum, '--m', 'Aluminum', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');
xline(glass, '--k', 'Glass', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');

% Change the color of each segment on the plot according to the material
hold on
for i = 1:length(distance_vec)-1
    if distance_vec(i+1) <= steel
        plot(distance_vec(i:i+1), average_temp(i:i+1), 'r')
    elseif distance_vec(i+1) <= granite
        plot(distance_vec(i:i+1), average_temp(i:i+1), 'b')
    elseif distance_vec(i+1) <= aluminum
        plot(distance_vec(i:i+1), average_temp(i:i+1), 'g')
    else
        plot(distance_vec(i:i+1), average_temp(i:i+1), 'm')
    end
end
hold off

% Question 2: Calculate the average slope of each material

% Calculate the average slope of each material
materials = {'Steel', 'Granite', 'Aluminum', 'Glass'};
boundaries = [0, steel, granite, aluminum, glass, max(distance_vec)];
average_slopes = zeros(1, length(materials));

for j = 1:length(materials)
    start_idx = find(distance_vec >= boundaries(j), 1);
    end_idx = find(distance_vec <= boundaries(j+1), 1, 'last');
    average_slopes(j) = (average_temp(end_idx) - average_temp(start_idx)) / (distance_vec(end_idx) - distance_vec(start_idx));
end

% Display the average slopes
for j = 1:length(materials)
    fprintf('Average slope for %s: %.2f\n', materials{j}, average_slopes(j));
end
% create the ratio of the average slopes against the steel slope
ratio = average_slopes / average_slopes(1);
for j = 2:length(materials)
    fprintf('Ratio of %s slope to steel slope: %.2f\n', materials{j}, ratio(j));
end

%% Section 2: Heat transfer through fins

% Q6: Calculate the average temperature for each fin

% Read the data from the fins and ribs CSV file
fins = readtable('fins&ribs.csv');
frame2dis =fins.Var1(end)/400; % Distance between frames in mm

% Extract the distance and temperature data for each fin
distance_vec = fins.Var1 ./ frame2dis; % Distance from base in mm
fin1_temp = [fins.Li1, fins.Li2, fins.Li7];
fin2_temp = [fins.Li4, fins.Li3, fins.Li8];
fin3_temp = [fins.Li6, fins.Li5, fins.Li9];

% Calculate the average temperature for each fin
average_temp_fin1 = mean(fin1_temp, 2);
average_temp_fin2 = mean(fin2_temp, 2);
average_temp_fin3 = mean(fin3_temp, 2);

% Plot the temperature vs. distance for Fin 1
figure;
plot(distance_vec, average_temp_fin1, 'r');
xlabel('Distance from base [mm]');
ylabel('Temperature [C]');
title('Temperature vs Distance from Base (Fin 1)');
grid on;

% Plot the temperature vs. distance for Fin 2
figure;
plot(distance_vec, average_temp_fin2, 'b');
xlabel('Distance from base [mm]');
ylabel('Temperature [C]');
title('Temperature vs Distance from Base (Fin 2)');
grid on;


% Plot the temperature vs. distance for Fin 3
figure;
plot(distance_vec, average_temp_fin3, 'g');
xlabel('Distance from base [mm]');
ylabel('Temperature [C]');
title('Temperature vs Distance from Base (Fin 3)');
grid on;


% Q7 and Q8: Plot the theoretical temperature distribution for each fin
L  = [100 200 400];
h = [5 3 2];
k = 18;
Tb  = 300;
T0 = 25;
averageTemp = {average_temp_fin1, average_temp_fin2, average_temp_fin3};

for i = 1:length(L)
    Tres = [];
    Tres2 = [];
    Tres3 = [];
    x = 0:0.1:L(i);
    m = sqrt(h(i)/(k*L(i)));
    for j = 1:length(x)
        Tres(j) = T0 + (Tb - T0) * (cosh(m * (L(i) - x(j))) + (h(i) / (m * k)) * sinh(m * (L(i) - x(j)))) / (cosh(m * L(i)) + (h(i) / (m * k)) * sinh(m * L(i)));
        Tres2(j) = T0 + (Tb - T0) * cosh(m * (L(i) - x(j))) / cosh(m * L(i));
        Tres3(j) = T0 + (Tb - T0) * exp(-m * x(j));
    end
    figure;
    plot(x, Tres, 'r', 'DisplayName', 'Theoretical');
    hold on;
    plot(x, Tres2, 'g', 'DisplayName', 'Adiabatic');
    plot(x, Tres3, 'c', 'DisplayName', 'Infinite');
    % plot(distance_vec, averageTemp{i}, 'b', 'DisplayName', 'Experimental');
    title(['Theoretical vs Experimental Temperature Distribution for Fin ', num2str(i)]);
    xlabel('Distance from base [mm]');
    ylabel('Temperature [C]');
    legend('show');
    grid on;
end


%% Section 3: Heat transfer through radiation

F = [0.057173 0.06307 0.068606 0.075026 0.082249 0.090503 0.099981];
T = [46.3 47.3 49.7 53 56.3 64.8 71.7];

sigma = 5.67*10^(-8);

T_kelv= T+273.15;

for i = 1:7
    h(i)=(sigma*(F(i)*(300+273.15)^(4)+(1-F(i))*(25+273.15)^(4)-(T_kelv(i))^4))/(T_kelv(i)-(25+273.15));
end
figure;
hold on;
plot(T,h,'r-')
grid on;
title('h(T) vs T')
xlabel('T[°C]')
ylabel('h[W/m^2K]')