pkg load statistics
clc;
clear all;
close all;

# ============================================================
# Lab 1 - Part 3
# Poisson Counting Process
# ============================================================

lambda = 10;                 # events / sec
mean_interarrival = 1 / lambda;

# ------------------------------------------------------------
# TASK A
# Generate 100 interarrival times and draw N(t)
# Interarrival times of a Poisson process are exponential
# ------------------------------------------------------------
number_of_events = 100;

interarrival_times = exprnd(mean_interarrival, 1, number_of_events);
arrival_times = cumsum(interarrival_times);

figure(1);
stairs([0 arrival_times], 0:number_of_events, "linewidth", 1.2);
grid on;
title("Poisson counting process for lambda = 10 events/sec");
xlabel("Time [sec]");
ylabel("N(t)");

fprintf("For a Poisson process with rate lambda = %.2f events/sec,\n", lambda);
fprintf("the interarrival times follow an exponential distribution with mean 1/lambda = %.4f sec.\n\n", mean_interarrival);


# ------------------------------------------------------------
# TASK B
# Estimate the average number of events per unit time
# for different numbers of generated events
# ------------------------------------------------------------
events_list = [100, 200, 400, 700, 1000, 5000, 50000, 100000];

observation_times = zeros(1, columns(events_list));
lambda_estimates = zeros(1, columns(events_list));

for i = 1:columns(events_list)
  n_events = events_list(i);

  interarrival_times_i = exprnd(mean_interarrival, 1, n_events);
  arrival_times_i = cumsum(interarrival_times_i);

  observation_times(i) = arrival_times_i(end);
  lambda_estimates(i) = n_events / observation_times(i);
endfor

fprintf("Estimated rate lambda_hat for increasing number of events:\n");
fprintf("--------------------------------------------------------\n");

for i = 1:columns(events_list)
  fprintf("n = %6d   total_time = %12.6f sec   lambda_hat = %10.6f events/sec\n", ...
          events_list(i), observation_times(i), lambda_estimates(i));
endfor

fprintf("\n");

figure(2);
plot(events_list, lambda_estimates, "o-", "linewidth", 1.5);
hold on;
plot(events_list, lambda * ones(size(events_list)), "--", "linewidth", 1.5);
hold off;
grid on;
set(gca, "xscale", "log");
title("Estimated rate versus number of generated events");
xlabel("Number of generated events (log scale)");
ylabel("Estimated rate lambda_hat");
legend("Estimated lambda_hat", "True lambda = 10");


# ------------------------------------------------------------
# Theoretical question for two non-overlapping windows
# DeltaT1 = 5 sec, DeltaT2 = 10 sec
# For a Poisson process:
# N(DeltaT) ~ Poisson(lambda * DeltaT)
# ------------------------------------------------------------
deltaT1 = 5;
deltaT2 = 10;

mean_events_window_1 = lambda * deltaT1;
mean_events_window_2 = lambda * deltaT2;
mean_total_events = mean_events_window_1 + mean_events_window_2;

fprintf("Theoretical results for two non-overlapping time windows:\n");
fprintf("--------------------------------------------------------\n");
fprintf("Window 1: DeltaT1 = %d sec  -> mean number of events = %.6f\n", ...
        deltaT1, mean_events_window_1);
fprintf("Window 2: DeltaT2 = %d sec  -> mean number of events = %.6f\n", ...
        deltaT2, mean_events_window_2);
fprintf("Total over both windows      -> mean number of events = %.6f\n", ...
        mean_total_events);

fprintf("\nDistribution remarks:\n");
fprintf("N1 ~ Poisson(lambda * DeltaT1) = Poisson(%.2f)\n", mean_events_window_1);
fprintf("N2 ~ Poisson(lambda * DeltaT2) = Poisson(%.2f)\n", mean_events_window_2);
fprintf("N1 + N2 ~ Poisson(%.2f)\n", mean_total_events);
