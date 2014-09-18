function EulerVerletComparrison( theta0, v0, T, dt )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Converts angle to radians
theta0 = theta0*(pi/180);

time = 0:dt:T;
l = length(time);
L = 1;	  % Length of pendulum in m
g = 9.8;  % m/s^2

v = zeros(1,l);
a = zeros(1,l);

%***********************************************************
% Euler Method

euler = zeros(1,l);

% First step of Euler's Method
euler(1) = theta0;
v(1) = v0;

for n = 1:l-1

euler(n+1) = euler(n) + v(n)*dt;
a(n) = (-g/L)*sin(euler(n));
v(n+1) = v(n) + a(n)*dt;

end

%***********************************************************
% Verlet Method

verlet = zeros(1,l);

verlet(1) = theta0;
a(1) = (-g/L)*sin(verlet(1));

% This estimates a theta before releasing the pendulum to start the Verlet method
earlyTheta = verlet(1) -v0*dt + .5*a(1)*dt^2;

verlet(2) = 2*verlet(1) - earlyTheta + dt^2*a(1);

for n = 2:l-1

a(n) = (-g/L)*sin(verlet(n));
verlet(n+1) = 2*verlet(n) - verlet(n-1) + dt^2*a(n);

end

%***********************************************************
% Comparrison

% For LaTeX use
%set(gca,'defaulttextinterpreter','latex');

% Change back into degrees

euler = euler.*(180/pi);
verlet = verlet.*(180/pi);

hold on

plot(time,euler,time,verlet,'linewidth',1)
legend('Euler','Verlet')
xlabel('time (s)')
ylabel('theta (degrees)')
title('Euler vs. Verlet Comparrison')

end