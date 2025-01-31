km = 0.0242; % Nm/A
Ra = 2.15; % Ohm
J = 0.0047; % kgm^2
f = 0.004; % Nms/rad
ng = 4.85; % gear ratio
s = tf('s');
P = km*ng/(Ra*J*s^2+Ra*f*s+ng^2*km^2*s);