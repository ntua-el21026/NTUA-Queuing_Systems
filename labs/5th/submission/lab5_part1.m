function lab5_part1()
% =========================================================================
% Lab 5 - Part A / Theme A
% Discrete Event Simulation of a Web Server
%
% File name:
%   lab5_part1.m
%
% Run in Octave with:
%   lab5_part1
%
% Implements:
%   A.1  Single DES run for lambda = 100 req/sec, plots and metrics.
%   A.2  Queue/system size plots for lambda in {50,100,150,200,250}.
%   A.3  95% confidence intervals from n >= 30 independent replications,
%        before and after transient removal.
%   A.4  Little's Law verification, with and without transient removal.
%
% Units:
%   - lambda is in requests/sec.
%   - service times and simulation time are in milliseconds.
%   - inter-arrival mean time is 1000/lambda ms.
%
% Model:
%   - Single FIFO server.
%   - Each request first becomes a type-1 job.
%   - After type-1 completion, the same request becomes a type-2 job and
%     re-enters the back of the same FIFO queue.
%   - After type-2 completion, the request leaves the system.
% =========================================================================

clc;
close all;

% Try to load the statistics package if available.
% The code does not require exprnd/lognrnd because local samplers are used.
% tinv/norminv are used if available; otherwise fallback approximations are used.
try
    pkg load statistics;
catch
    printf('Note: statistics package was not loaded. Local samplers/fallback quantiles will be used where needed.\n');
end_try_catch

% Fresh RNG state at the start of the script.
rand('state', sum(100 * clock));
randn('state', sum(100 * clock) + 12345);

% =========================================================================
% User-adjustable simulation parameters
% =========================================================================

max_req_A1          = 1e4;
max_req_A2          = 1e4;
n_runs_A3           = 30;
n_req_rep_A3        = 1e4;
warmup_fraction_A3  = 0.30;     % MQS-guided conservative transient removal
warmup_fraction_A4  = 0.30;     % same rule for Little's Law comparison

% =========================================================================
% A.1 - Single run for lambda = 100 requests/sec
% =========================================================================

printf('\n============================================================\n');
printf('A.1 - Single DES run for lambda = 100 req/sec\n');
printf('============================================================\n');

lambda_A1 = 100;
r = sim_webserver(lambda_A1, max_req_A1);

R1_A1 = mean(r.resp_t1);
R2_A1 = mean(r.resp_t2);
R_A1  = mean(r.total_resp);

tp1_A1 = length(r.resp_t1) / (r.t_end / 1000);
tp2_A1 = length(r.resp_t2) / (r.t_end / 1000);

N1_A1 = tw_mean(r.log_q1, r.log_t);
N2_A1 = tw_mean(r.log_q2, r.log_t);
N_A1  = tw_mean(r.log_q1 + r.log_q2, r.log_t);

printf('Simulated %d fully served requests up to t = %.4f ms = %.4f sec\n', ...
       r.n_served, r.t_end, r.t_end / 1000);
printf('Total arrivals observed:        %d\n', r.n_arrived);
printf('R1, type-1 phase response:      %.6f ms\n', R1_A1);
printf('R2, type-2 phase response:      %.6f ms\n', R2_A1);
printf('R end-to-end:                   %.6f ms\n', R_A1);
printf('Type-1 jobs served per second:  %.6f\n', tp1_A1);
printf('Type-2 jobs served per second:  %.6f\n', tp2_A1);
printf('Mean #type-1 in system, tw:     %.6f\n', N1_A1);
printf('Mean #type-2 in system, tw:     %.6f\n', N2_A1);
printf('Mean total #jobs in system, tw: %.6f\n', N_A1);

% A.1 Plot 1: cumulative arrivals and fully served requests.
figure(1);
clf;
plot(r.log_t / 1000, r.log_arrived, 'LineWidth', 1.2);
hold on;
plot(r.log_t / 1000, r.log_served, 'LineWidth', 1.2);
grid on;
xlabel('Time (sec)');
ylabel('Cumulative requests');
title('A.1 - Arrived and fully served requests, lambda = 100 req/sec');
legend('Arrived requests', 'Fully served requests', 'Location', 'northwest');

% A.1 Plot 2: number of type-1 and type-2 jobs in the system.
figure(2);
clf;

subplot(2, 1, 1);
plot(r.log_t / 1000, r.log_q1, 'LineWidth', 1.0);
grid on;
xlabel('Time (sec)');
ylabel('Type-1 jobs');
title('A.1 - Type-1 jobs in system, queue + service');

subplot(2, 1, 2);
plot(r.log_t / 1000, r.log_q2, 'LineWidth', 1.0);
grid on;
xlabel('Time (sec)');
ylabel('Type-2 jobs');
title('A.1 - Type-2 jobs in system, queue + service');

% A.1 Plot 3: response times and throughput as bar plots.
figure(3);
clf;

subplot(1, 2, 1);
bar([R1_A1, R2_A1]);
grid on;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Type 1', 'Type 2'});
ylabel('Mean response time (ms)');
title('A.1 - Mean response time per phase');

subplot(1, 2, 2);
bar([tp1_A1, tp2_A1]);
grid on;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Type 1', 'Type 2'});
ylabel('Completed jobs per second');
title('A.1 - Mean served jobs per second');

% =========================================================================
% A.2 - Stationarity across lambda values
% =========================================================================

printf('\n============================================================\n');
printf('A.2 - Stationarity across lambda values\n');
printf('============================================================\n');

S1_mean = exp(1 + 0.5);        % E[S1] = exp(mu + sigma^2/2), mu = sigma = 1
S2_mean = 0.75;                % E[S2] = (0.5 + 1) / 2
S_total = S1_mean + S2_mean;

lambda_stability_bound = 1000 / S_total;

printf('E[S1] = exp(1.5) = %.6f ms\n', S1_mean);
printf('E[S2] = %.6f ms\n', S2_mean);
printf('E[S1] + E[S2] = %.6f ms\n', S_total);
printf('Analytical stability condition: rho = lambda*(E[S1]+E[S2])/1000 < 1\n');
printf('Stability bound: lambda < %.6f req/sec\n', lambda_stability_bound);

lambdas_A2 = [50, 100, 150, 200, 250];
runs_A2 = cell(1, length(lambdas_A2));

figure(4);
clf;

for k = 1:length(lambdas_A2)
    lam = lambdas_A2(k);
    rho = lam * S_total / 1000;

    if rho < 1
        status = 'stable';
    else
        status = 'unstable';
    endif

    printf('Running lambda = %d req/sec, rho = %.6f -> %s ... ', lam, rho, status);
    fflush(stdout);

    runs_A2{k} = sim_webserver(lam, max_req_A2);

    printf('done, t_end = %.4f sec, arrivals = %d, served = %d\n', ...
           runs_A2{k}.t_end / 1000, runs_A2{k}.n_arrived, runs_A2{k}.n_served);

    q_total = runs_A2{k}.log_q1 + runs_A2{k}.log_q2;

    subplot(length(lambdas_A2), 1, k);
    plot(runs_A2{k}.log_t / 1000, q_total, 'LineWidth', 1.0);
    grid on;
    xlabel('Time (sec)');
    ylabel('Total jobs');
    title(sprintf('A.2 - lambda = %d req/sec, rho = %.4f, %s', lam, rho, status));
endfor

printf('\nAnalytical classification for A.2:\n');
for k = 1:length(lambdas_A2)
    lam = lambdas_A2(k);
    rho = lam * S_total / 1000;

    if rho < 1
        printf('  lambda = %3d req/sec: rho = %.6f < 1 -> stationary\n', lam, rho);
    else
        printf('  lambda = %3d req/sec: rho = %.6f >= 1 -> non-stationary\n', lam, rho);
    endif
endfor

% =========================================================================
% A.3 - Confidence intervals before and after transient removal
% =========================================================================

printf('\n============================================================\n');
printf('A.3 - Confidence intervals from independent replications\n');
printf('============================================================\n');

lambdas_A3 = [100, 180];

CI1_mean_all = zeros(length(lambdas_A3), 2);
CI1_med_all  = zeros(length(lambdas_A3), 2);
CI2_mean_all = zeros(length(lambdas_A3), 2);
CI2_med_all  = zeros(length(lambdas_A3), 2);

CI1_mean_tr = zeros(length(lambdas_A3), 2);
CI1_med_tr  = zeros(length(lambdas_A3), 2);
CI2_mean_tr = zeros(length(lambdas_A3), 2);
CI2_med_tr  = zeros(length(lambdas_A3), 2);

% Pilot MQS plots for cutoff visualization.
figure(5);
clf;

for li = 1:length(lambdas_A3)
    lam = lambdas_A3(li);

    pilot = sim_webserver(lam, n_req_rep_A3);
    q_pilot_total = pilot.log_q1 + pilot.log_q2;
    mqs_pilot = running_tw_mean(q_pilot_total, pilot.log_t);
    cutoff_pilot = warmup_fraction_A3 * pilot.t_end;

    subplot(length(lambdas_A3), 1, li);
    plot(pilot.log_t / 1000, mqs_pilot, 'LineWidth', 1.0);
    hold on;
    yl = ylim();
    plot([cutoff_pilot / 1000, cutoff_pilot / 1000], yl, '--', 'LineWidth', 1.0);
    grid on;
    xlabel('Time (sec)');
    ylabel('MQS(t)');
    title(sprintf('A.3 - MQS(t), lambda = %d req/sec, cutoff = %.0f%% of run time', ...
                  lam, 100 * warmup_fraction_A3));
    legend('MQS(t)', 'Transient cutoff', 'Location', 'southeast');
endfor

for li = 1:length(lambdas_A3)
    lam = lambdas_A3(li);

    N1_all = zeros(n_runs_A3, 1);
    N2_all = zeros(n_runs_A3, 1);
    N1_tr  = zeros(n_runs_A3, 1);
    N2_tr  = zeros(n_runs_A3, 1);

    printf('\nlambda = %d req/sec: running %d independent replications ', ...
           lam, n_runs_A3);

    for rep = 1:n_runs_A3
        rk = sim_webserver(lam, n_req_rep_A3);

        N1_all(rep) = tw_mean(rk.log_q1, rk.log_t);
        N2_all(rep) = tw_mean(rk.log_q2, rk.log_t);

        cutoff_ms = warmup_fraction_A3 * rk.t_end;

        N1_tr(rep) = tw_mean_interval(rk.log_q1, rk.log_t, cutoff_ms, rk.t_end);
        N2_tr(rep) = tw_mean_interval(rk.log_q2, rk.log_t, cutoff_ms, rk.t_end);

        printf('.');
        fflush(stdout);
    endfor

    printf(' done\n');

    [ci1_mean_all, ci1_med_all] = calc_ci(N1_all, 0.05);
    [ci2_mean_all, ci2_med_all] = calc_ci(N2_all, 0.05);

    [ci1_mean_tr, ci1_med_tr] = calc_ci(N1_tr, 0.05);
    [ci2_mean_tr, ci2_med_tr] = calc_ci(N2_tr, 0.05);

    CI1_mean_all(li, :) = ci1_mean_all;
    CI1_med_all(li, :)  = ci1_med_all;
    CI2_mean_all(li, :) = ci2_mean_all;
    CI2_med_all(li, :)  = ci2_med_all;

    CI1_mean_tr(li, :) = ci1_mean_tr;
    CI1_med_tr(li, :)  = ci1_med_tr;
    CI2_mean_tr(li, :) = ci2_mean_tr;
    CI2_med_tr(li, :)  = ci2_med_tr;

    printf('  Full trace, type 1: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           ci1_mean_all(1), ci1_mean_all(2), ci1_med_all(1), ci1_med_all(2));
    printf('  Full trace, type 2: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           ci2_mean_all(1), ci2_mean_all(2), ci2_med_all(1), ci2_med_all(2));

    printf('  After transient removal, type 1: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           ci1_mean_tr(1), ci1_mean_tr(2), ci1_med_tr(1), ci1_med_tr(2));
    printf('  After transient removal, type 2: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           ci2_mean_tr(1), ci2_mean_tr(2), ci2_med_tr(1), ci2_med_tr(2));
endfor

% =========================================================================
% A.4 - Little's Law verification for lambda = 100 requests/sec
% =========================================================================

printf('\n============================================================\n');
printf('A.4 - Little''s Law verification for lambda = 100 req/sec\n');
printf('============================================================\n');

cutoff_A4 = warmup_fraction_A4 * r.t_end;

[N_sim_full, R_full, N_little_full, err_full] = little_stats(r, lambda_A1, 0);
[N_sim_tr, R_tr, N_little_tr, err_tr] = little_stats(r, lambda_A1, cutoff_A4);

printf('Without transient removal:\n');
printf('  R = %.6f ms\n', R_full);
printf('  N_sim, time-weighted total jobs = %.6f\n', N_sim_full);
printf('  lambda * R = (%.6f req/ms) * %.6f ms = %.6f\n', ...
       lambda_A1 / 1000, R_full, N_little_full);
printf('  Percent error = %.6f %%\n', err_full);

printf('\nWith transient removal:\n');
printf('  Cutoff = %.6f ms = %.6f sec, %.0f%% of run time\n', ...
       cutoff_A4, cutoff_A4 / 1000, 100 * warmup_fraction_A4);
printf('  R = %.6f ms\n', R_tr);
printf('  N_sim, time-weighted total jobs after cutoff = %.6f\n', N_sim_tr);
printf('  lambda * R = (%.6f req/ms) * %.6f ms = %.6f\n', ...
       lambda_A1 / 1000, R_tr, N_little_tr);
printf('  Percent error = %.6f %%\n', err_tr);

% =========================================================================
% Final summary
% =========================================================================

printf('\n============================================================\n');
printf('FINAL SUMMARY - Lab 5, Part A\n');
printf('============================================================\n');

printf('\nA.1  lambda = %d req/sec, %d fully served requests\n', lambda_A1, max_req_A1);
printf('  t_end = %.6f ms = %.6f sec\n', r.t_end, r.t_end / 1000);
printf('  arrivals = %d, fully served = %d\n', r.n_arrived, r.n_served);
printf('  R1 = %.6f ms, R2 = %.6f ms, R = %.6f ms\n', R1_A1, R2_A1, R_A1);
printf('  throughput: type-1 = %.6f jobs/sec, type-2 = %.6f jobs/sec\n', tp1_A1, tp2_A1);
printf('  time-weighted mean jobs: N1 = %.6f, N2 = %.6f, N = %.6f\n', ...
       N1_A1, N2_A1, N_A1);

printf('\nA.2  Stability bound\n');
printf('  lambda < %.6f req/sec\n', lambda_stability_bound);

for k = 1:length(lambdas_A2)
    lam = lambdas_A2(k);
    rho = lam * S_total / 1000;

    if rho < 1
        printf('  lambda = %3d: rho = %.6f -> stationary\n', lam, rho);
    else
        printf('  lambda = %3d: rho = %.6f -> non-stationary\n', lam, rho);
    endif
endfor

printf('\nA.3  95%% confidence intervals, n = %d replications\n', n_runs_A3);
printf('  Transient removal rule: remove first %.0f%% of each run, guided by MQS(t) plots.\n', ...
       100 * warmup_fraction_A3);

for li = 1:length(lambdas_A3)
    lam = lambdas_A3(li);

    printf('\n  lambda = %d req/sec\n', lam);

    printf('    Full trace:\n');
    printf('      type 1: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           CI1_mean_all(li, 1), CI1_mean_all(li, 2), ...
           CI1_med_all(li, 1),  CI1_med_all(li, 2));
    printf('      type 2: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           CI2_mean_all(li, 1), CI2_mean_all(li, 2), ...
           CI2_med_all(li, 1),  CI2_med_all(li, 2));

    printf('    After transient removal:\n');
    printf('      type 1: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           CI1_mean_tr(li, 1), CI1_mean_tr(li, 2), ...
           CI1_med_tr(li, 1),  CI1_med_tr(li, 2));
    printf('      type 2: mean CI [%.6f, %.6f], median CI [%.6f, %.6f]\n', ...
           CI2_mean_tr(li, 1), CI2_mean_tr(li, 2), ...
           CI2_med_tr(li, 1),  CI2_med_tr(li, 2));
endfor

printf('\nA.4  Little''s Law, lambda = %d req/sec\n', lambda_A1);
printf('  Without transient removal:\n');
printf('    R = %.6f ms, N_sim = %.6f, lambda*R = %.6f, error = %.6f %%\n', ...
       R_full, N_sim_full, N_little_full, err_full);
printf('  With transient removal:\n');
printf('    R = %.6f ms, N_sim = %.6f, lambda*R = %.6f, error = %.6f %%\n', ...
       R_tr, N_sim_tr, N_little_tr, err_tr);

printf('\nDone.');

endfunction

% =========================================================================
% Local functions
% =========================================================================

function res = sim_webserver(lambda, max_req)
    % Run one DES simulation of the web server.
    %
    % lambda  : external request arrival rate in requests/sec
    % max_req : stop after this many requests have fully finished
    %
    % Returned fields:
    %   log_t, log_q1, log_q2, log_arrived, log_served
    %   resp_t1, resp_t2, total_resp
    %   completed_depart2, completed_arrival
    %   n_arrived, n_served, t_end

    ARRIVAL = 1;
    DEPART1 = 2;
    DEPART2 = 3;

    inter_arr_mean = 1000 / lambda;       % ms

    % Event list rows: [time, event_type]
    events = [sample_exp_mean(inter_arr_mean), ARRIVAL];

    % FIFO queue rows: [job_type, phase_start_time, request_id]
    queue = [];

    busy = false;
    cur_job = [0, 0, 0];                  % [job_type, phase_start_time, request_id]

    % Logs start at t = 0 with an empty system.
    log_t       = 0;
    log_q1      = 0;
    log_q2      = 0;
    log_arrived = 0;
    log_served  = 0;

    resp_t1 = [];
    resp_t2 = [];
    total_resp = [];

    arrival_times = [];
    depart1_times = [];

    completed_ids = [];
    completed_arrival = [];
    completed_depart2 = [];

    n_arrived = 0;
    n_served  = 0;
    t = 0;

    while n_served < max_req
        if isempty(events)
            break;
        endif

        % Pick the event with the smallest timestamp.
        [t, idx] = min(events(:, 1));
        etype = events(idx, 2);
        events(idx, :) = [];

        if etype == ARRIVAL
            n_arrived = n_arrived + 1;
            req_id = n_arrived;

            arrival_times(req_id) = t;

            % Schedule the next external arrival.
            events = [events; t + sample_exp_mean(inter_arr_mean), ARRIVAL];

            if ~busy
                busy = true;
                cur_job = [1, t, req_id];
                events = [events; t + sample_lognormal(1, 1), DEPART1];
            else
                queue = [queue; 1, t, req_id];
            endif

        elseif etype == DEPART1
            req_id = cur_job(3);

            depart1_times(req_id) = t;
            resp_t1 = [resp_t1, t - cur_job(2)];

            % The same request becomes a type-2 job at the back of the FIFO queue.
            queue = [queue; 2, t, req_id];

            % Start the next job from the FIFO queue.
            [cur_job, queue] = pop_next(queue);

            if cur_job(1) == 1
                events = [events; t + sample_lognormal(1, 1), DEPART1];
            else
                events = [events; t + sample_uniform(0.5, 1.0), DEPART2];
            endif

        elseif etype == DEPART2
            req_id = cur_job(3);

            resp_t2 = [resp_t2, t - cur_job(2)];
            total_resp = [total_resp, t - arrival_times(req_id)];

            completed_ids = [completed_ids, req_id];
            completed_arrival = [completed_arrival, arrival_times(req_id)];
            completed_depart2 = [completed_depart2, t];

            n_served = n_served + 1;

            if ~isempty(queue)
                [cur_job, queue] = pop_next(queue);

                if cur_job(1) == 1
                    events = [events; t + sample_lognormal(1, 1), DEPART1];
                else
                    events = [events; t + sample_uniform(0.5, 1.0), DEPART2];
                endif
            else
                busy = false;
                cur_job = [0, 0, 0];
            endif
        endif

        % Count type-1 and type-2 jobs in the system: FIFO queue + current service.
        if isempty(queue)
            nq1 = 0;
            nq2 = 0;
        else
            nq1 = sum(queue(:, 1) == 1);
            nq2 = sum(queue(:, 1) == 2);
        endif

        nq1 = nq1 + (cur_job(1) == 1);
        nq2 = nq2 + (cur_job(1) == 2);

        log_t       = [log_t,       t];
        log_q1      = [log_q1,      nq1];
        log_q2      = [log_q2,      nq2];
        log_arrived = [log_arrived, n_arrived];
        log_served  = [log_served,  n_served];
    endwhile

    res.log_t       = log_t;
    res.log_q1      = log_q1;
    res.log_q2      = log_q2;
    res.log_arrived = log_arrived;
    res.log_served  = log_served;

    res.resp_t1 = resp_t1;
    res.resp_t2 = resp_t2;
    res.total_resp = total_resp;

    res.completed_ids = completed_ids;
    res.completed_arrival = completed_arrival;
    res.completed_depart2 = completed_depart2;

    res.n_arrived = n_arrived;
    res.n_served  = n_served;
    res.t_end = t;
endfunction

function [job, queue_out] = pop_next(queue)
    % FIFO pop.
    job = queue(1, :);
    queue_out = queue(2:end, :);
endfunction

function x = sample_exp_mean(mean_value)
    % Exponential random variable with the given mean.
    u = rand();

    while u <= 0
        u = rand();
    endwhile

    x = -mean_value * log(u);
endfunction

function x = sample_lognormal(mu, sigma)
    % LogNormal(mu, sigma), where mu/sigma are the parameters of the
    % underlying normal distribution.
    x = exp(mu + sigma * randn());
endfunction

function x = sample_uniform(a, b)
    % Uniform(a,b).
    x = a + (b - a) * rand();
endfunction

function avg = tw_mean(vals, times)
    % Time-weighted average of a step function.
    % vals(i) is assumed to hold over [times(i), times(i+1)).
    avg = tw_mean_interval(vals, times, times(1), times(end));
endfunction

function avg = tw_mean_interval(vals, times, t_start, t_stop)
    % Time-weighted average of a step function over [t_start, t_stop].
    vals = vals(:)';
    times = times(:)';

    if length(times) < 2
        avg = NaN;
        return;
    endif

    t_start = max(t_start, times(1));
    t_stop  = min(t_stop,  times(end));

    if t_stop <= t_start
        avg = NaN;
        return;
    endif

    area = 0;
    total_time = 0;

    for i = 1:(length(times) - 1)
        a = max(times(i),     t_start);
        b = min(times(i + 1), t_stop);

        if b > a
            area = area + vals(i) * (b - a);
            total_time = total_time + (b - a);
        endif
    endfor

    if total_time <= 0
        avg = NaN;
    else
        avg = area / total_time;
    endif
endfunction

function mqs = running_tw_mean(vals, times)
    % Running time-weighted average:
    % MQS(t_j) = (1 / (t_j - t_0)) * integral_{t_0}^{t_j} Q(s) ds.
    vals = vals(:)';
    times = times(:)';

    n = length(times);
    mqs = zeros(1, n);

    if n < 2
        return;
    endif

    area = 0;
    t0 = times(1);
    mqs(1) = vals(1);

    for i = 1:(n - 1)
        dt = times(i + 1) - times(i);

        if dt < 0
            error('Times must be non-decreasing.');
        endif

        area = area + vals(i) * dt;
        denom = times(i + 1) - t0;

        if denom > 0
            mqs(i + 1) = area / denom;
        else
            mqs(i + 1) = mqs(i);
        endif
    endfor
endfunction

function [ci_mean, ci_med] = calc_ci(xs, alpha)
    % Confidence intervals for mean and median.
    %
    % Mean CI:
    %   xbar +/- t_{n-1,1-alpha/2} * s/sqrt(n)
    %
    % Median CI:
    %   order-statistics interval [X_(k), X_(n+1-k)] with
    %   k = floor(n/2 - z*sqrt(n)/2).

    if nargin < 2
        alpha = 0.05;
    endif

    xs = xs(:);
    xs = xs(~isnan(xs));
    xs = sort(xs);

    n = length(xs);

    if n < 2
        ci_mean = [NaN, NaN];
        ci_med  = [NaN, NaN];
        return;
    endif

    xbar = mean(xs);
    s = std(xs);

    eta = t_quantile(1 - alpha / 2, n - 1);
    ci_mean = [xbar - eta * s / sqrt(n), xbar + eta * s / sqrt(n)];

    z = normal_quantile(1 - alpha / 2);
    k = floor(n / 2 - z * sqrt(n) / 2);
    k = max(k, 1);
    k = min(k, floor((n + 1) / 2));

    ci_med = [xs(k), xs(n + 1 - k)];
endfunction

function q = normal_quantile(p)
    % Normal quantile. Uses norminv if available, otherwise erfinv formula.
    if exist('norminv', 'file') || exist('norminv', 'builtin')
        q = norminv(p);
    else
        q = sqrt(2) * erfinv(2 * p - 1);
    endif
endfunction

function q = t_quantile(p, df)
    % Student-t quantile.
    % Uses tinv if available. If not available, uses standard normal as a
    % fallback, except for df = 29 where the common n = 30 value is included.
    if exist('tinv', 'file') || exist('tinv', 'builtin')
        q = tinv(p, df);
    else
        if df == 29 && abs(p - 0.975) < 1e-12
            q = 2.045230;
        else
            q = normal_quantile(p);
        endif
    endif
endfunction

function [N_sim, R_mean, N_little, err_percent] = little_stats(r, lambda, cutoff_ms)
    % Little's Law verification for one trace.
    %
    % N_sim:
    %   time-weighted mean total number of jobs in the system.
    %
    % R_mean:
    %   mean end-to-end response time of requests completed after cutoff.
    %
    % N_little:
    %   lambda_ms * R_mean.
    %
    % err_percent:
    %   relative absolute error between N_little and N_sim.

    q_total = r.log_q1 + r.log_q2;

    if cutoff_ms <= 0
        N_sim = tw_mean(q_total, r.log_t);
        idx = true(size(r.total_resp));
    else
        N_sim = tw_mean_interval(q_total, r.log_t, cutoff_ms, r.t_end);
        idx = (r.completed_arrival >= cutoff_ms) & (r.completed_depart2 >= cutoff_ms);
    endif

    if sum(idx) == 0
        R_mean = NaN;
        N_little = NaN;
        err_percent = NaN;
        return;
    endif

    R_mean = mean(r.total_resp(idx));

    lambda_ms = lambda / 1000;
    N_little = lambda_ms * R_mean;

    if abs(N_sim) < eps
        err_percent = NaN;
    else
        err_percent = 100 * abs(N_little - N_sim) / abs(N_sim);
    endif
endfunction
