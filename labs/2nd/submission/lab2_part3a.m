clc;
clear all;
close all;
format long g;

pkg load queueing;

% Part 3: Birth-death process for M/M/1/3
% lambda_i = lambda / (i + 2), i = 0,1,2
% mu_i = mu, independent of state
% lambda = 5 customers/sec
% mu = 10 customers/sec
% Initial state: empty system

lambda = 5;
mu = 10;

states = [0, 1, 2, 3];
K = 3;

% Initial probability vector: system starts empty
initial_state = [1, 0, 0, 0];

% Birth and death rates
% State mapping:
% 0 -> 1 : lambda/2
% 1 -> 2 : lambda/3
% 2 -> 3 : lambda/4
births_B = [lambda / 2, lambda / 3, lambda / 4];
deaths_D = [mu, mu, mu];

% -------------------------------------------------------------------------
% (i) Transition-rate matrix (generator matrix)
% -------------------------------------------------------------------------
transition_matrix = ctmcbd(births_B, deaths_D);

fprintf("Transition-rate matrix Q:\n");
disp(transition_matrix);

% -------------------------------------------------------------------------
% (ii) Ergodic probabilities with ctmc(Q)
%     Also compute the theoretical probabilities from balance equations
% -------------------------------------------------------------------------
P_octave = ctmc(transition_matrix);

ratios = zeros(1, K + 1);
ratios(1) = 1;

for n = 2 : (K + 1)
  ratios(n) = ratios(n - 1) * births_B(n - 1) / deaths_D(n - 1);
endfor

P_theory = ratios / sum(ratios);

fprintf("Theoretical stationary probabilities:\n");
disp(P_theory);

fprintf("Octave stationary probabilities:\n");
disp(P_octave);

fprintf("Comparison table [state, P_theory, P_octave, abs_error]:\n");
comparison_table = [states(:), P_theory(:), P_octave(:), abs(P_theory(:) - P_octave(:))];
disp(comparison_table);

% Plot stationary probabilities
figure(1);
bar(states, P_octave, 0.5);
grid on;
xlabel("State");
ylabel("Probability");
title("M/M/1/3 - Stationary Probabilities");

% -------------------------------------------------------------------------
% (iii) Blocking probability
% -------------------------------------------------------------------------
blocking_probability = P_octave(end);

fprintf("Blocking probability P_block = P3 = %.12f\n", blocking_probability);

% -------------------------------------------------------------------------
% (iv) Mean number of customers in the system
% -------------------------------------------------------------------------
mean_clients = sum(states .* P_octave);

fprintf("Mean number of customers E[n(t)] = %.12f\n", mean_clients);

% -------------------------------------------------------------------------
% (v) Mean number of customers served in 60 sec after equilibrium
% Throughput at equilibrium = sum over states of (service rate in state i)*P_i
% For states 1,2,3 service rate is mu; for state 0 it is 0
% -------------------------------------------------------------------------
service_rates_by_state = [0, deaths_D];
throughput_equilibrium = sum(service_rates_by_state .* P_octave);
avg_served_60 = throughput_equilibrium * 60;

fprintf("Equilibrium throughput = %.12f customers/sec\n", throughput_equilibrium);
fprintf("Mean customers served in 60 sec = %.12f\n", avg_served_60);

% -------------------------------------------------------------------------
% (vi) Transient probabilities until all state probabilities are within 1%
%      of the stationary probabilities
% -------------------------------------------------------------------------
tolerance = 0.01;
dt = 0.01;
tmax = 20;

time_values = 0 : dt : tmax;
n_times = numel(time_values);
n_states = numel(states);

prob_history = zeros(n_times, n_states);
convergence_index = n_times;

for idx = 1 : n_times
  t = time_values(idx);
  P_t = ctmc(transition_matrix, t, initial_state);
  prob_history(idx, :) = P_t;

  if all(abs(P_t - P_octave) < tolerance)
    convergence_index = idx;
    break;
  endif
endfor

time_values = time_values(1 : convergence_index);
prob_history = prob_history(1 : convergence_index, :);

convergence_time = time_values(end);

fprintf("Time until all probabilities are within 1%% of steady state: %.4f sec\n", convergence_time);

figure(2);

for s = 1 : n_states
  subplot(2, 2, s);
  plot(time_values, prob_history(:, s), "linewidth", 1.5);
  hold on;
  plot(time_values, P_octave(s) * ones(size(time_values)), "--", "linewidth", 1.0);
  grid on;
  xlabel("Time (sec)");
  ylabel(sprintf("P_%d(t)", states(s)));
  title(sprintf("State %d", states(s)));
endfor
