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
