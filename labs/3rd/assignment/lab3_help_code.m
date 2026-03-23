clear; clc;
rand("seed", 1);

lambda = 5;
mu = 5;
K = 10;
T = 10000;
B = 1000;

state = 0;

arrival_attempts = 0;
rejections = 0;
sum_states = 0;
state_counts = zeros(1, K+1);

arrival_attempts_B = 0;
rejections_B = 0;
sum_states_B = 0;
state_counts_B = zeros(1, K+1);

for step = 1:T

  % count current state
  state_counts(state + 1) = state_counts(state + 1) + 1;
  sum_states = sum_states + state;

  if step > B
    state_counts_B(state + 1) = state_counts_B(state + 1) + 1;
    sum_states_B = sum_states_B + state;
  endif

  % choose event
  if state == 0
    event = "arrival";
    arrival_attempts = arrival_attempts + 1;
    state_new = 1;

    if step > B
      arrival_attempts_B = arrival_attempts_B + 1;
    endif

  else
    if rand() < lambda / (lambda + mu)
      % arrival attempt
      arrival_attempts = arrival_attempts + 1;
      if step > B
        arrival_attempts_B = arrival_attempts_B + 1;
      endif

      if state < K
        event = "arrival";
        state_new = state + 1;
      else
        event = "rejection";
        state_new = state;
        rejections = rejections + 1;
        if step > B
          rejections_B = rejections_B + 1;
        endif
      endif

    else
      event = "departure";
      state_new = state - 1;
    endif
  endif

  if step <= 30
    printf("Step %2d | state=%2d | event=%s | new state=%2d\n", ...
           step, state, event, state_new);
  endif

  state = state_new;
endfor

P = state_counts / sum(state_counts);
L = sum_states / T;
P_rej = rejections / arrival_attempts;

P_B = state_counts_B / sum(state_counts_B);
L_B = sum_states_B / (T - B);
P_rej_B = rejections_B / arrival_attempts_B;

fprintf('\nUsing all transitions:\n');
fprintf('Mean number in system = %f\n', L);
fprintf('Rejection probability = %f\n', P_rej);
disp('State probabilities:');
disp(P);

fprintf('\nIgnoring first B transitions:\n');
fprintf('Mean number in system = %f\n', L_B);
fprintf('Rejection probability = %f\n', P_rej_B);
disp('State probabilities after B:');
disp(P_B);
