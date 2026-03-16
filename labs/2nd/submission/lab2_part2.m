clc;
clear all;
close all;
format long g;

pkg load queueing;

% Part 2: M/M/1 analysis with Octave
% Arrival rate: lambda = 10 customers/min
% Service rate: mu in (10, 20] customers/min

lambda = 10;                     % customers/min
mu_min = lambda + 0.05;          % avoid the singular point mu = lambda
mu_max = 20;
mu_step = 0.05;
mu_values = mu_min : mu_step : mu_max;

n_points = numel(mu_values);

U  = zeros(1, n_points);   % utilization
R  = zeros(1, n_points);   % average response time E(T)
Qm = zeros(1, n_points);   % average number of customers in the system
X  = zeros(1, n_points);   % throughput
p0 = zeros(1, n_points);   % probability system is empty

for i = 1 : n_points
  mu = mu_values(i);
  [U(i), R(i), Qm(i), X(i), p0(i)] = qsmm1(lambda, mu);
endfor

fprintf("Accepted service rates for ergodicity:\n");
fprintf("mu > lambda = %.2f customers/min, with mu <= %.2f customers/min\n\n", lambda, mu_max);
fprintf("So the acceptable range is: %.2f < mu <= %.2f customers/min\n\n", lambda, mu_max);

% Figure 1: Utilization vs service rate
figure(1);
plot(mu_values, U, "linewidth", 1.6);
grid on;
xlabel("Service rate mu (customers/min)");
ylabel("Utilization U");
title("M/M/1 - Utilization vs Service Rate");

% Figure 2: Average response time vs service rate
figure(2);
plot(mu_values, R, "linewidth", 1.6);
grid on;
xlabel("Service rate mu (customers/min)");
ylabel("Average response time E(T) (min)");
title("M/M/1 - Average Response Time vs Service Rate");

% Figure 3: Average number of customers vs service rate
figure(3);
plot(mu_values, Qm, "linewidth", 1.6);
grid on;
xlabel("Service rate mu (customers/min)");
ylabel("Average number of customers E[n(t)]");
title("M/M/1 - Average Number of Customers vs Service Rate");

% Figure 4: Throughput vs service rate
figure(4);
plot(mu_values, X, "linewidth", 1.6);
grid on;
xlabel("Service rate mu (customers/min)");
ylabel("Throughput X (customers/min)");
title("M/M/1 - Throughput vs Service Rate");

% Optional terminal outputs for quick inspection
[~, idx_min_delay] = min(R);

fprintf("Minimum delay in the scanned range occurs at:\n");
fprintf("mu = %.2f customers/min\n", mu_values(idx_min_delay));
fprintf("E(T) = %.6f min\n", R(idx_min_delay));
fprintf("Throughput there is X = %.6f customers/min\n", X(idx_min_delay));
