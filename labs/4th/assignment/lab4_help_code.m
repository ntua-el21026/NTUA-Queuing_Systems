pkg load queueing
lambda = 100; % Calls per hour
mu = 12;      % Calls per agent per hour
m = 10;       % Number of agents

A = lambda / mu; % Offered Load (8.33 Erlangs)
Pw = erlangc(A, m); % Probability of waiting
