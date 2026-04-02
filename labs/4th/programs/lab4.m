clear;
clc;
close all;
pkg load queueing
format long g

% Lab 4 - Erlang C
% Minimal extension of the provided helper code

lambda = 100;           % Calls per hour
h_minutes = 3;          % Mean call duration [minutes]
h = h_minutes / 60;     % Mean call duration [hours]
mu = 1 / h;             % Calls per server per hour = 20
m = 10;                 % Number of agents / servers

% Service Level threshold t
% The statement leaves it as "t".
% Default choice here: 20 sec (typical 80/20 target from the intro slides).
% If your instructor wants another value, change only this line.
t_seconds = 20;
t = t_seconds / 3600;   % [hours]

% a) Offered load
A = lambda * h;         % Equivalent to lambda / mu

% Erlang C stability condition
if A >= m
  error('Unstable Erlang-C system: A must be smaller than m.');
endif

% b) P0
sum_terms = 0;
for k = 0:(m - 1)
  sum_terms = sum_terms + (A^k) / factorial(k);
endfor
tail_term = ((A^m) / factorial(m)) * (m / (m - A));
P0 = 1 / (sum_terms + tail_term);

% c) Probability of waiting
Pw = erlangc(A, m);          % From queueing package
Pw_formula = tail_term * P0; % Closed-form cross-check

% d) Service Level
% SL(t) = 1 - Pw * exp(-(m*mu - lambda) * t)
SL = 1 - Pw * exp(-(m * mu - lambda) * t);

% e) Average Speed of Answer
ASA_hours = (Pw * h) / (m - A);
ASA_minutes = ASA_hours * 60;
ASA_seconds = ASA_hours * 3600;

% Extra sanity metric
rho = A / m;  % Occupancy / utilization

fprintf('================ Erlang C - Lab 4 ================\n');
fprintf('lambda                  = %.10f calls/hour\n', lambda);
fprintf('h                       = %.10f minutes/call\n', h_minutes);
fprintf('mu                      = %.10f calls/server/hour\n', mu);
fprintf('m                       = %d servers\n', m);
fprintf('t                       = %.10f seconds\n', t_seconds);
fprintf('---------------------------------------------------\n');
fprintf('(a) Offered load A      = %.10f Erlangs\n', A);
fprintf('(b) P0                  = %.10f\n', P0);
fprintf('(c) Pw (erlangc)        = %.10f\n', Pw);
fprintf('    Pw (closed-form)    = %.10f\n', Pw_formula);
fprintf('    abs difference      = %.10e\n', abs(Pw - Pw_formula));
fprintf('(d) Service Level SL(t) = %.10f\n', SL);
fprintf('                        = %.10f %%\n', 100 * SL);
fprintf('(e) ASA                 = %.10f hours\n', ASA_hours);
fprintf('                        = %.10f minutes\n', ASA_minutes);
fprintf('                        = %.10f seconds\n', ASA_seconds);
fprintf('Extra: rho = A / m      = %.10f\n', rho);
fprintf('                        = %.10f %%\n', 100 * rho);
fprintf('===================================================\n');
