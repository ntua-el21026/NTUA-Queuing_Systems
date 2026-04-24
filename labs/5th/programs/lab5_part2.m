function lab5_part2()
% =========================================================================
% Lab 5 - Part B
% Open Jackson Queueing Network
%
% File name:
%   lab5_part2.m
%
% Run in Octave with:
%   lab5_part2
%
% Implements:
%   B.2  intensities(lambda1, lambda2, mu1, ..., mu5)
%   B.3  mean_clients(lambda1, lambda2, mu1, ..., mu5)
%   B.4  numerical evaluation for the given parameters
%   B.5  bottleneck queue and maximum lambda1 for ergodicity
%   B.6  end-to-end mean delay plot versus lambda1
%
% Network routing used:
%   Q1 -> Q2     with probability 1/2
%   Q1 -> Q3     with probability 1/4
%   Q1 -> exit   with probability 1/4
%
%   Q2 -> Q3     with probability 1/3
%   Q2 -> Q4     with probability 2/3
%
%   Q3 -> Q5     with probability 1
%
%   Q4 -> Q5     with probability 3/5
%   Q4 -> exit   with probability 2/5
%
%   Q5 -> exit   with probability 1
%
% Effective arrival rates:
%   Gamma1 = lambda1
%   Gamma2 = lambda2 + (1/2)lambda1
%   Gamma3 = (5/12)lambda1 + (1/3)lambda2
%   Gamma4 = (1/3)lambda1 + (2/3)lambda2
%   Gamma5 = (37/60)lambda1 + (11/15)lambda2
%
% All rates are in clients/sec, so delays are in seconds.
% =========================================================================

clc;
close all;
format long g;

% -------------------------------------------------------------------------
% Given numerical parameters from the statement
% -------------------------------------------------------------------------

lambda1 = 10;
lambda2 = 6;

mu1 = 15;
mu2 = 20;
mu3 = 10;
mu4 = 14;
mu5 = 12;

mu = [mu1, mu2, mu3, mu4, mu5];

fprintf('\n============================================================\n');
fprintf('Lab 5 - Part B: Open Jackson Queueing Network\n');
fprintf('============================================================\n');

fprintf('\nInput parameters:\n');
fprintf('  lambda1 = %.10f clients/sec\n', lambda1);
fprintf('  lambda2 = %.10f clients/sec\n', lambda2);
fprintf('  mu1     = %.10f clients/sec\n', mu1);
fprintf('  mu2     = %.10f clients/sec\n', mu2);
fprintf('  mu3     = %.10f clients/sec\n', mu3);
fprintf('  mu4     = %.10f clients/sec\n', mu4);
fprintf('  mu5     = %.10f clients/sec\n', mu5);

% -------------------------------------------------------------------------
% B.2 - Effective arrival rates and intensities
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('B.2 - Effective arrival rates and traffic intensities\n');
fprintf('============================================================\n');

gamma = effective_arrival_rates(lambda1, lambda2);

fprintf('\nEffective arrival rates Gamma_i:\n');
for i = 1:5
    fprintf('  Gamma%d = %.10f clients/sec\n', i, gamma(i));
endfor

[rho, ergodic] = intensities(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5);

fprintf('\nErgodicity flag: %d\n', ergodic);
if ergodic == 1
    fprintf('The network is ergodic because all rho_i are smaller than 1.\n');
else
    fprintf('The network is not ergodic because at least one rho_i is >= 1.\n');
endif

% -------------------------------------------------------------------------
% B.3 - Mean number of clients in each queue
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('B.3 - Mean number of clients in each queue\n');
fprintf('============================================================\n');

N = mean_clients(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5);

fprintf('\nMean clients E[N_i] in each queue:\n');
for i = 1:5
    fprintf('  E[N%d] = %.10f\n', i, N(i));
endfor

N_total = sum(N);

fprintf('\nTotal mean number of clients in the network:\n');
fprintf('  E[N] = sum_i E[N_i] = %.10f\n', N_total);

% -------------------------------------------------------------------------
% B.4 - End-to-end mean delay for the given parameters
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('B.4 - End-to-end mean delay for the given parameters\n');
fprintf('============================================================\n');

external_rate = lambda1 + lambda2;
T_end_to_end = end_to_end_delay(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5);

fprintf('\nTotal external arrival rate:\n');
fprintf('  gamma_external = lambda1 + lambda2 = %.10f clients/sec\n', external_rate);

fprintf('\nEnd-to-end mean delay by Little''s Law:\n');
fprintf('  E[T] = E[N] / gamma_external = %.10f seconds\n', T_end_to_end);

% -------------------------------------------------------------------------
% B.5 - Bottleneck and maximum lambda1 for ergodicity
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('B.5 - Bottleneck queue and maximum lambda1\n');
fprintf('============================================================\n');

[max_rho, bottleneck_current] = max(rho);

fprintf('\nCurrent bottleneck queue, based on maximum rho_i:\n');
fprintf('  bottleneck = Q%d\n', bottleneck_current);
fprintf('  max rho    = %.10f\n', max_rho);

[lambda1_max, limiting_queue, bounds] = lambda1_max_for_ergodicity(lambda2, mu1, mu2, mu3, mu4, mu5);

fprintf('\nUpper bounds on lambda1 imposed by each queue:\n');
for i = 1:5
    fprintf('  From Q%d: lambda1 < %.10f\n', i, bounds(i));
endfor

fprintf('\nMaximum lambda1 for ergodicity with lambda2 fixed at %.10f:\n', lambda2);
fprintf('  lambda1_max = %.10f clients/sec\n', lambda1_max);
fprintf('  limiting queue = Q%d\n', limiting_queue);

if lambda1 < lambda1_max
    fprintf('  Current lambda1 = %.10f is below lambda1_max, so the current network is ergodic.\n', lambda1);
else
    fprintf('  Current lambda1 = %.10f is not below lambda1_max, so the current network is not ergodic.\n', lambda1);
endif

% -------------------------------------------------------------------------
% B.6 - Plot end-to-end delay versus lambda1
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('B.6 - End-to-end delay plot versus lambda1\n');
fprintf('============================================================\n');

lambda1_values = linspace(0.1 * lambda1_max, 0.99 * lambda1_max, 300);
delay_values = zeros(size(lambda1_values));

for k = 1:length(lambda1_values)
    delay_values(k) = end_to_end_delay(lambda1_values(k), lambda2, ...
                                       mu1, mu2, mu3, mu4, mu5);
endfor

fprintf('\nPlot range:\n');
fprintf('  lambda1_min = %.10f clients/sec\n', lambda1_values(1));
fprintf('  lambda1_max_plot = %.10f clients/sec\n', lambda1_values(end));
fprintf('  lambda1_max_ergodic = %.10f clients/sec\n', lambda1_max);

fprintf('\nDelay values at selected points:\n');

selected_factors = [0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99];
for k = 1:length(selected_factors)
    l1_sel = selected_factors(k) * lambda1_max;
    d_sel = end_to_end_delay(l1_sel, lambda2, mu1, mu2, mu3, mu4, mu5);

    fprintf('  lambda1 = %.2f * lambda1_max = %.10f -> E[T] = %.10f sec\n', ...
            selected_factors(k), l1_sel, d_sel);
endfor

figure(1);
clf;

plot(lambda1_values, delay_values, 'LineWidth', 1.5);
grid on;
xlabel('lambda_1 (clients/sec)');
ylabel('End-to-end mean delay E[T] (sec)');
title('B.6 - End-to-end mean delay versus lambda_1');

hold on;
yl = ylim();
plot([lambda1, lambda1], yl, '--', 'LineWidth', 1.0);
plot([lambda1_max, lambda1_max], yl, ':', 'LineWidth', 1.0);
legend('E[T]', 'Given lambda_1 = 10', 'lambda_{1,max}', 'Location', 'northwest');
hold off;

try
    print(gcf, 'lab5_part2_delay.png', '-dpng', '-r200');
    fprintf('\nSaved figure as lab5_part2_delay.png\n');
catch
    fprintf('\nCould not save figure automatically. Please export it manually if needed.\n');
end_try_catch

% -------------------------------------------------------------------------
% Final summary for easy copy into report
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('FINAL SUMMARY - Lab 5, Part B\n');
fprintf('============================================================\n');

fprintf('\nRouting equations:\n');
fprintf('  Gamma1 = lambda1\n');
fprintf('  Gamma2 = lambda2 + (1/2)lambda1\n');
fprintf('  Gamma3 = (5/12)lambda1 + (1/3)lambda2\n');
fprintf('  Gamma4 = (1/3)lambda1 + (2/3)lambda2\n');
fprintf('  Gamma5 = (37/60)lambda1 + (11/15)lambda2\n');

fprintf('\nGiven parameters:\n');
fprintf('  lambda1 = %.10f, lambda2 = %.10f\n', lambda1, lambda2);
fprintf('  mu = [%.10f %.10f %.10f %.10f %.10f]\n', mu1, mu2, mu3, mu4, mu5);

fprintf('\nEffective arrival rates:\n');
fprintf('  Gamma = [');
fprintf(' %.10f', gamma);
fprintf(' ]\n');

fprintf('\nTraffic intensities:\n');
fprintf('  rho = [');
fprintf(' %.10f', rho);
fprintf(' ]\n');
fprintf('  ergodic = %d\n', ergodic);

fprintf('\nMean clients:\n');
fprintf('  E[N_i] = [');
fprintf(' %.10f', N);
fprintf(' ]\n');
fprintf('  E[N] = %.10f\n', N_total);

fprintf('\nEnd-to-end delay:\n');
fprintf('  E[T] = %.10f seconds\n', T_end_to_end);

fprintf('\nBottleneck:\n');
fprintf('  current bottleneck = Q%d, rho = %.10f\n', bottleneck_current, max_rho);
fprintf('  lambda1_max = %.10f clients/sec\n', lambda1_max);
fprintf('  limiting queue for lambda1_max = Q%d\n', limiting_queue);

fprintf('\nDone.');

endfunction

% =========================================================================
% Required function for B.2
% =========================================================================

function [rho, ergodic] = intensities(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5, verbose)
    % Computes traffic intensities rho_i for the five queues.
    %
    % Required public use:
    %   [rho, ergodic] = intensities(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5)
    %
    % Optional internal use:
    %   verbose = false suppresses printing.

    if nargin < 8
        verbose = true;
    endif

    mu = [mu1, mu2, mu3, mu4, mu5];
    gamma = effective_arrival_rates(lambda1, lambda2);

    rho = gamma ./ mu;
    ergodic = all(rho < 1);

    if verbose
        fprintf('\nTraffic intensities rho_i = Gamma_i / mu_i:\n');
        for i = 1:5
            fprintf('  rho%d = %.10f\n', i, rho(i));
        endfor
    endif
endfunction

% =========================================================================
% Required function for B.3
% =========================================================================

function N = mean_clients(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5)
    % Computes the mean number of clients in each M/M/1 queue:
    %   E[N_i] = rho_i / (1 - rho_i)
    %
    % If the system is not ergodic, returns NaN values.

    [rho, ergodic] = intensities(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5, false);

    if ergodic == 0
        N = NaN(1, 5);
        return;
    endif

    N = rho ./ (1 - rho);
endfunction

% =========================================================================
% Helper functions
% =========================================================================

function gamma = effective_arrival_rates(lambda1, lambda2)
    % Effective arrival rates at each queue, derived from the routing graph.

    gamma1 = lambda1;
    gamma2 = lambda2 + (1/2) * lambda1;
    gamma3 = (5/12) * lambda1 + (1/3) * lambda2;
    gamma4 = (1/3) * lambda1 + (2/3) * lambda2;
    gamma5 = (37/60) * lambda1 + (11/15) * lambda2;

    gamma = [gamma1, gamma2, gamma3, gamma4, gamma5];
endfunction

function T = end_to_end_delay(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5)
    % Computes end-to-end mean delay using Little's Law:
    %   E[T] = E[N] / (lambda1 + lambda2)

    N = mean_clients(lambda1, lambda2, mu1, mu2, mu3, mu4, mu5);
    external_rate = lambda1 + lambda2;

    if any(isnan(N)) || external_rate <= 0
        T = NaN;
        return;
    endif

    T = sum(N) / external_rate;
endfunction

function [lambda1_max, limiting_queue, bounds] = lambda1_max_for_ergodicity(lambda2, mu1, mu2, mu3, mu4, mu5)
    % Computes the maximum lambda1 such that all queues remain ergodic,
    % with lambda2 and all mu_i fixed.
    %
    % Conditions:
    %   Gamma1 < mu1
    %   Gamma2 < mu2
    %   Gamma3 < mu3
    %   Gamma4 < mu4
    %   Gamma5 < mu5
    %
    % where:
    %   Gamma1 = lambda1
    %   Gamma2 = lambda2 + (1/2)lambda1
    %   Gamma3 = (5/12)lambda1 + (1/3)lambda2
    %   Gamma4 = (1/3)lambda1 + (2/3)lambda2
    %   Gamma5 = (37/60)lambda1 + (11/15)lambda2

    bound1 = mu1;
    bound2 = 2 * (mu2 - lambda2);
    bound3 = (12/5) * (mu3 - (1/3) * lambda2);
    bound4 = 3 * (mu4 - (2/3) * lambda2);
    bound5 = (60/37) * (mu5 - (11/15) * lambda2);

    bounds = [bound1, bound2, bound3, bound4, bound5];

    [lambda1_max, limiting_queue] = min(bounds);

    if lambda1_max < 0
        lambda1_max = NaN;
    endif
endfunction
