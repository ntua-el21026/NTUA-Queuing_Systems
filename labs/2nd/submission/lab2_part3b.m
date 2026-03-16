clc;
clear all;
close all;
format long g;

pkg load queueing;

% Part 3 - subquestion (vii)
% Repeat the 60-sec serviced-customers calculation for:
% (i) lambda = 5, mu = 1
% (ii) lambda = 5, mu = 5
% (iii) lambda = 5, mu = 20
%
% Also inspect how stationary probabilities and convergence speed change.

lambda = 5;
mu_values = [1, 5, 20];

states = [0, 1, 2, 3];
K = 3;
initial_state = [1, 0, 0, 0];

tolerance = 0.01;
dt = 0.01;
tmax = 20;

summary_results = zeros(numel(mu_values), 6);
% columns:
% [mu, P0, P1, P2, P3, avg_served_60]

for case_idx = 1 : numel(mu_values)

  mu = mu_values(case_idx);

  births_B = [lambda / 2, lambda / 3, lambda / 4];
  deaths_D = [mu, mu, mu];

  transition_matrix = ctmcbd(births_B, deaths_D);
  P_stationary = ctmc(transition_matrix);

  % Throughput and mean customers served in 60 sec
  service_rates_by_state = [0, deaths_D];
  throughput_equilibrium = sum(service_rates_by_state .* P_stationary);
  avg_served_60 = throughput_equilibrium * 60;

  % Transient analysis until all state probabilities are within 1%
  time_values = 0 : dt : tmax;
  n_times = numel(time_values);
  n_states = numel(states);

  prob_history = zeros(n_times, n_states);
  convergence_index = n_times;

  for idx = 1 : n_times
    t = time_values(idx);
    P_t = ctmc(transition_matrix, t, initial_state);
    prob_history(idx, :) = P_t;

    if all(abs(P_t - P_stationary) < tolerance)
      convergence_index = idx;
      break;
    endif
  endfor

  time_values = time_values(1 : convergence_index);
  prob_history = prob_history(1 : convergence_index, :);
  convergence_time = time_values(end);

  fprintf("\n============================================================\n");
  fprintf("Case %d: lambda = %.2f, mu = %.2f\n", case_idx, lambda, mu);
  fprintf("Stationary probabilities [P0 P1 P2 P3]:\n");
  disp(P_stationary);
  fprintf("Blocking probability = %.12f\n", P_stationary(end));
  fprintf("Equilibrium throughput = %.12f customers/sec\n", throughput_equilibrium);
  fprintf("Mean customers served in 60 sec = %.12f\n", avg_served_60);
  fprintf("Convergence time (within 1%% of steady state) = %.4f sec\n", convergence_time);

  summary_results(case_idx, :) = [mu, P_stationary, avg_served_60];

  % Plot probabilities for this case
  figure(case_idx);

  for s = 1 : n_states
    subplot(2, 2, s);
    plot(time_values, prob_history(:, s), "linewidth", 1.5);
    hold on;
    plot(time_values, P_stationary(s) * ones(size(time_values)), "--", "linewidth", 1.0);
    grid on;
    xlabel("Time (sec)");
    ylabel(sprintf("P_%d(t)", states(s)));
    title(sprintf("State %d, mu = %.2f", states(s), mu));
  endfor

endfor

fprintf("\n============================================================\n");
fprintf("Summary table [mu, P0, P1, P2, P3, avg_served_60]\n");
disp(summary_results);
