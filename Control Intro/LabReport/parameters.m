km = 0.0243; % Nm/A
Ra = 2.08; % Ohm
J = 0.0047; % kgm^2
f = 0.0077* (pi/180); % Nms/Deg
ng = 5.25; % gear ratio
s = tf('s');
P = km*ng/s*(Ra*J*s+Ra*f+ng^2*km^2);