km = 0.0243; % Nm/A
Ra = 2.08; % Ohm
J = 0.0047; % kgm^2
f = 0.0077; % Nms/sec
ng = 5.25; % gear ratio

Kp = 0.0485; %P controller
Ki = 3; %I controller

s = tf('s');
P = (ng*km*(180/pi)/(s*(Ra*J*s + (Ra*f + ng^2*km^2))));