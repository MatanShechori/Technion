%% Known Var of the Question

m = 1; % kg
g = 9.81; %m/sec^2
L0 = 0.2; %m
a = 1; % m
k = 4; %N/m
c1 = 0.1;
c3 = 0.001;
teta = 30*pi/180; % rad

%% quesion 5

f = @(x) k*(1- L0/(sqrt(x^2+a^2)))*x-m*g*sin(teta);

root  = fzero(f,1);
fprintf('root' , root);


%% question 6
L0 = 5; %m new initial length
sol = fzero(f,1);

disp(sol);

%% question 7



%% Q8


%% Q9


%% Q10

%% Q11


%% Q12

