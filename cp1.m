clc;        % Clear command line
clear;      % Clear workspace variables
close all   % Closes open figures 
%% Part 2

%{
    Define the f(t, y) function for Part 2.
    Input:
    t: Time step, a scalar
    y: Numerical solution of the current step, 2-element vector, [displacement, velocity]
    Output:
    A 2-element vector as derivative of y
%}

%--------------------------------------------------------------------------
% TODO: Implement the function f(t,y)
%--------------------------------------------------------------------------

w1 = 2;
z1 = 0.05;
k1 = 2;
f=@(t,y)[y(2), sin(k1*w1*t)-(2*z1*w1*y(2))-(w1^2)*y(1)];
N = 81;
y = zeros(N,2);
y(1,:) = [1,-1];
times = linspace(0.0, 8.0, N);
h = 0.1;
for i = 1:N-1
t = i*h;
k1 = f(t, y(i,:));
k2 = f(t+h/2, y(i,:)+k1*h/2);
y(i+1,:) = y(i,:) + k2*h;
disp( y(i,:))
end
disp( y(i+1,:))
%---------------------------------------------------------------------
% Set up for Part 2
%---------------------------------------------------------------------

tf = 8;                    % Final time for part 2
y0 = [0.,0.];              % Initial condition
%--------------------------------------------------------------------------
% TODO: Compute the analytical solution
% HINT: 1. Use the vector calculation wisely to avoid for-loops
%       2. Define variables for coefficients to simplify the expressions
%--------------------------------------------------------------------------
ya1 = 0;    % The displacement
ya2 = 0;    % The velocity
x_0 = 0;
xdot1_0 = 0;
w1 = 2;
z1 = 0.05;
k1 = 2;
wd = w1*(sqrt(1-(z1^2)));
% Constants
format long
H = (w1^2)-((k1^2)*(w1^2));
E = 2*z1*(w1^2)*(k1);
A1 = -E/((H^2)+(E^2));
B1 = H/((H^2)+(E^2));
C1 = x_0-A1;
C2 = ((xdot1_0) + (z1 * w1 * (x_0 - A1)) - (B1 * k1 * w1))/(wd);


% Analytic Solution
N = 81;
anal_arr = zeros(N,2);
anal_arr(1,:) = [1,-1];
for i2= 1:N-1
t2 = 0.1*i2;
x = exp(-z1*w1*t2)*(C1*cos(wd*t2) + C2*sin(wd*t2)) + A1*cos(k1*w1*t2) +B1*sin(k1*w1*t2);
anal_arr(i2+1,1) = x;
x_prime = (-z1)*w1*exp(-z1*w1*t2)*(C1*cos(wd*t2)+C2*sin(wd*t2)) +wd*exp(-z1*w1*t2)*(-C1*sin(wd*t2)+C2*cos(wd*t2)) + k1*w1*B1*cos(k1*w1*t2) -k1*w1*A1*sin(k1*w1*t2);
anal_arr(i2+1,2) = x_prime;
end

%--------------------------------------------------------------------------
% Part 2(1) & (2)
%--------------------------------------------------------------------------
% This generates the numerical solution using RK2
% function [t, y] = rk2(y0,f1 tf, dt);
f=@(t,y)[y(2), sin(2*w1*t)-(w1^2)*y(1)-(2*z1*w1*y(2))];
% f=@(t,y)[y(2), -(w1^2)*y(1)-(2*z1*w1*y(2))];
y = zeros(N,2);
y(1,:) = [1,-1];
times = linspace(0.0, 8.0, N);
h = 0.1;
for i = 1:N-1
t = i*h;
k1 = f(t, y(i,:));
k2 = f(t+h/2, y(i,:)+k1*h/2);
y(i+1,:) = y(i,:) + k2*h;
end


% Generate the numerical solution using ode45 in MATLAB.
t3 = 0;
y3 = 0;

%--------------------------------------------------------------------------
% TODO: Process and plot the results
%--------------------------------------------------------------------------
disp( y(i+1,:))
disp(anal_arr(N,:))
plot(times, y(:,1))
hold on
plot(times, anal_arr(:,1))
hold on
plot(times, y(:,2))
hold on
plot(times, anal_arr(:,2))

%--------------------------------------------------------------------------
% Part 2(3)
%--------------------------------------------------------------------------                                  
Nh = 8;                            % The number of step sizes to check
yr = [ ya1(end); ya2(end) ];       % Exact solution at last time step

hs = zeros(Nh,1);                  % Allocate the array for step sizes
es = zeros(Nh,1);                  % Allocate the array for errors.
for ii = 1:Nh
   plt = 0.1/2^(ii-1);            % Halving the step size each time 

    %----------------------------------------------------------------------
    % TODO: Implement the rest of the loop
    %----------------------------------------------------------------------
    % You need to call rk2 using the step size h
    % And then record the step size (in hs) and the numerical error (in es)
end

%--------------------------------------------------------------------------
% TODO: Process and plot the results
%--------------------------------------------------------------------------
% Here is an example code for the log-log plot revealing the order of
% accuracy
% er = es(end)*(hs/hs(end)).^2;
% fig2 = figure(2);
% plt = loglog(hs, es, 'b-',...
%        hs, er, 'r--');
% legend([plt(1), plt(2)], {'RK2', '2nd-order reference'},...
%         'Location', 'SouthEast')
% xlabel('Step Size')
% ylabel('Error')
% exportgraphics(fig2, './pics/p23.png', 'Resolution', 300,...
%                                        'BackgroundColor', 'white');
%% Part 3 (1)

z3 = 0.05;  % Initial zeta for part 3
o3 = 100;   % Omega for part 3
h = 0.001;  % Time step for part 3   

f3 = @(t,y) [0.0; 0.0];   % Function to be integrated for part 3


%--------------------------------------------------------------------------
% Setup for Part 3
%--------------------------------------------------------------------------

tf3 = 0.2;      % Final time for Part 3
y0 = [0.,0.];   % Time step for Part 3
T0 = 1;         % T0 for part 3
zmax = 2.0;     % Max zeta value
dz = 0.1;       % Time steps for zeta values

%--------------------------------------------------------------------------
% TODO: Analyze, design and plot the results.
%--------------------------------------------------------------------------

f=@(t,y)[y(2), sin(k*w1*t)-(2*z1*w1*y(2))-(w1^2)*y(1)];
N = 81;
y = zeros(N,2);
y(1,:) = [1,-1];
anal_arr = zeros(N,2);
anal_arr(1,:) = [1,-1];
t = linspace(0.0, 8.0, N);
error_arr = [0.2,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001; 0,0,0,0,0,0,0,0];
for count = 1:8
h = error_arr(1,count);
for i = 1:N-1
t = i*h;
k1 = f(t, y(i,:));
k2 = f(t+h/2, y(i,:)+k1*h/2);
y(i+1,:) = y(i,:) + k2*h;
end
for i2= 1:N-1
t2 = h*i2;
x = exp(-z1*w1*t2)*(C1*cos(wd*t2) + C2*sin(wd*t2)) + A1*cos(k*w1*t2) + B1*sin(k*w1*t2);
anal_arr(i2+1,1) = x;
x_prime = (-z1)*w1*exp(-z1*w1*t2)*(C1*cos(wd*t2)+C2*sin(wd*t2)) +wd*exp(-z1*w1*t2)*(-C1*sin(wd*t2)+C2*cos(wd*t2)) + k*w1*B1*cos(k*w1*t2) -k*w1*A1*sin(k*w1*t2);
anal_arr(i2+1,2) = x_prime;
end
error_total = abs(((anal_arr(i2+1,1))^2 + (anal_arr(i2+1,2))^2)^0.5 - ((y(i+1,1))^2 +(y(i+1,2))^2)^0.5);
error_arr(2,count) = (error_total);
end
disp(error_arr)
%% Problem 3_1
w1 = 2;
z1 = 0;
k1 = 2;
f=@(t,y)[y(2), -(2*z1*w1*y(2))-(w1^2)*y(1)];
N = 401;
y = zeros(N,2);
y(1,:) = [1,-1];
times = linspace(0.0, 40.0, N);
h = 0.1;
% RK2
for i = 1:N-1
t = i*h;
k1 = f(t, y(i,:));
k2 = f(t+h/2, y(i,:)+k1*h/2);
y(i+1,:) = y(i,:) + k2*h;
yt = y(i,:);



end

plot(times, y(:,1))
hold on
plot(times, y(:,2))



for i = 1:N-1
t = i*h;
k1 = f(t, y(i,:));
k2 = f(t+(h/2), y(i,:)+(k1*h/2));
k3= f(t+(h/2), y(i,:)+(k2*h/2));
k4= f(t+h, y(i,:)+(k3*h));
y(i+1,:) = y(i,:) + ((1/6)*(k1+(2*k2)+(2*k3)+k4)*h);
rt = y(i,:);
% disp("RK4 y(i,:) = ")
% disp( y(i,:))
end
% disp("RK4 y(i+1,:) = ")
% disp( y(i+1,:))
plot(times, y(:,1))
hold on
plot(times, y(:,2))
