% Lab 5 - Theme A: Web Server DES
% Demo.
%
% Setup: one FIFO queue, one server. Each request becomes a type-1 job;
% when it finishes, the same request becomes a type-2 job at the back of
% the queue (not the front). When the type-2 finishes, the request leaves.
%
% What you still need to do:
%   A.1  just the plots. The numerical results are computed and printed
%        by the script.
%   A.2  fill in the lambda values in lambdas_A2 and add the plots.
%   A.3  plot MQS(t), pick a cutoff, and redo the CIs after removing
%        the transient. See the TODO inside the loop.
%   A.4  verify lambda * R = N, with R = mean(resp_t1) + mean(resp_t2).
%
% Response times: each request has two phases. The type-1 phase starts at
% arrival; the type-2 phase starts later, when the type-1 finishes. So
% resp_t1 and resp_t2 are the two phase durations separately, and the
% end-to-end R is their sum -- that is the R for Little's Law.
%
% Units: lambda in req/sec, service times in ms, so inter-arrival is
% 1000/lambda ms. Keep everything in ms inside the sim.
%
% Useful functions from the statistics package:
%   exprnd(mean)       exponential samples, for inter-arrival times
%   lognrnd(mu, sigma) lognormal. Note: mu and sigma are the parameters of
%                      the underlying normal, so E[X] = exp(mu + sigma^2/2),
%                      not exp(mu).
%   tinv(p, df)        t-distribution quantile (mean CI)
%   norminv(p)         normal quantile (median CI)

clear; clc;
pkg load statistics;

% fresh RNG state at the start of each run
rand('state', sum(100*clock));
randn('state', sum(100*clock));

function res = sim_webserver(lambda, max_req)
    % Run one simulation.
    % lambda  : arrivals per second
    % max_req : stop after this many requests have fully finished
    %
    % Returns a struct with:
    %   log_t        timestamps of each event (ms)
    %   log_q1       type-1 jobs in the system at each event
    %   log_q2       type-2 jobs in the system at each event
    %   log_arrived  cumulative arrivals at each event
    %   log_served   cumulative fully-served requests at each event
    %   resp_t1      per type-1 job, time from phase start to DEPART1 (ms)
    %   resp_t2      per type-2 job, time from phase start to DEPART2 (ms)
    %   n_arrived    total arrivals
    %   n_served     total fully served requests
    %   t_end        time of the last event (ms)

    ARRIVAL = 1; DEPART1 = 2; DEPART2 = 3;
    inter_arr = 1000 / lambda;

    % event list, each row is [time, type]. Start with one arrival scheduled.
    events = [exprnd(inter_arr), ARRIVAL];

    % queue, each row is [job_type, phase_start_time]. We store the phase
    % start time so resp_ti is simply t_departure - phase_start for both
    % phases.
    queue = [];

    busy = false;
    cur_job = [0 0];   % [type, phase_start_time]

    log_t = []; log_q1 = []; log_q2 = [];
    log_arrived = []; log_served = [];
    resp_t1 = []; resp_t2 = [];
    n_arrived = 0;
    n_served  = 0;
    t = 0;

    while n_served < max_req
        if isempty(events), break; end

        % pick the event with the smallest time and remove it
        [t, idx] = min(events(:,1));
        etype = events(idx, 2);
        events(idx, :) = [];

        if etype == ARRIVAL
            n_arrived = n_arrived + 1;
            events = [events; t + exprnd(inter_arr), ARRIVAL];

            if ~busy
                % server free, start this request's type-1 now
                busy = true;
                cur_job = [1, t];
                events = [events; t + lognrnd(1, 1), DEPART1];
            else
                % server busy, put it in the queue as type-1
                queue = [queue; 1, t];
            end

        elseif etype == DEPART1
            % type-1 just finished. Record its time and create the type-2
            % at the back of the queue with phase start = now.
            resp_t1 = [resp_t1, t - cur_job(2)];
            queue = [queue; 2, t];

            % pop the head of the queue and start serving it
            [cur_job, queue] = pop_next(queue);
            if cur_job(1) == 1
                events = [events; t + lognrnd(1, 1), DEPART1];
            else
                events = [events; t + 0.5 + rand * 0.5, DEPART2];
            end

        elseif etype == DEPART2
            % type-2 just finished, the request leaves the system
            resp_t2 = [resp_t2, t - cur_job(2)];
            n_served = n_served + 1;

            if ~isempty(queue)
                [cur_job, queue] = pop_next(queue);
                if cur_job(1) == 1
                    events = [events; t + lognrnd(1, 1), DEPART1];
                else
                    events = [events; t + 0.5 + rand * 0.5, DEPART2];
                end
            else
                busy = false;
                cur_job = [0 0];
            end
        end

        % count type-1 and type-2 jobs in the system (queue + current)
        if isempty(queue)
            nq1 = 0; nq2 = 0;
        else
            nq1 = sum(queue(:,1) == 1);
            nq2 = sum(queue(:,1) == 2);
        end
        nq1 = nq1 + (cur_job(1) == 1);
        nq2 = nq2 + (cur_job(1) == 2);

        log_t       = [log_t,       t];
        log_q1      = [log_q1,      nq1];
        log_q2      = [log_q2,      nq2];
        log_arrived = [log_arrived, n_arrived];
        log_served  = [log_served,  n_served];
    end

    res.log_t       = log_t;
    res.log_q1      = log_q1;
    res.log_q2      = log_q2;
    res.log_arrived = log_arrived;
    res.log_served  = log_served;
    res.resp_t1     = resp_t1;
    res.resp_t2     = resp_t2;
    res.n_arrived   = n_arrived;
    res.n_served    = n_served;
    res.t_end       = t;
end

function [job, queue_out] = pop_next(queue)
    % take the first job off the queue (FIFO)
    job = queue(1, :);
    queue_out = queue(2:end, :);
end

function avg = tw_mean(vals, times)
    % Time-weighted average of a step function: vals(i) holds from
    % times(i) to times(i+1). Needed for the queue size, where a plain
    % mean() would over-weight the busy periods (events happen more
    % often when the queue is big).
    dt = diff(times);
    avg = sum(vals(1:end-1) .* dt) / sum(dt);
end

function [ci_mean, ci_med] = calc_ci(xs, alpha)
    % 95% CIs for mean and median of xs. Mean uses the t-distribution;
    % median uses order statistics.
    if nargin < 2, alpha = 0.05; end
    n  = length(xs);
    xs = sort(xs(:));

    xbar = mean(xs);
    s    = std(xs);
    eta  = tinv(1 - alpha/2, n - 1);
    ci_mean = [xbar - eta*s/sqrt(n), xbar + eta*s/sqrt(n)];

    z = norminv(1 - alpha/2);
    k = max(floor(n/2 - z*sqrt(n)/2), 1);
    ci_med = [xs(k), xs(n+1-k)];
end

% A.1: single run at lambda = 100.
% The script runs the simulation and computes all the numerical results.
% You only need to add the plots using the data stored in r: r.log_t,
% r.log_arrived, r.log_served, r.log_q1, r.log_q2.

lambda  = 100;
max_req = 1e4;
r = sim_webserver(lambda, max_req);

printf('Simulated %d requests up to t = %.1f ms\n', r.n_served, r.t_end);
printf('R1 (type-1 phase):             %.4f ms\n', mean(r.resp_t1));
printf('R2 (type-2 phase):             %.4f ms\n', mean(r.resp_t2));
printf('R end-to-end:                  %.4f ms\n', mean(r.resp_t1) + mean(r.resp_t2));

% mean number of type-i jobs served per second
tp1 = length(r.resp_t1) / (r.t_end / 1000);
tp2 = length(r.resp_t2) / (r.t_end / 1000);
printf('Type-1 jobs served per second: %.4f\n', tp1);
printf('Type-2 jobs served per second: %.4f\n', tp2);

% time-weighted mean number of type-i jobs in the system
N1_single = tw_mean(r.log_q1, r.log_t);
N2_single = tw_mean(r.log_q2, r.log_t);
printf('Mean #type-1 in system (tw):   %.4f\n', N1_single);
printf('Mean #type-2 in system (tw):   %.4f\n', N2_single);

% A.2: stationarity across different lambda values.
% Fill in lambdas_A2 below with the values you want to study, then plot
% the queue size over time for each. The analytical bound is printed
% for comparison.

% Analytical bound: each request occupies the server for E[S1]+E[S2] ms,
% so the system is stationary as long as lambda*(E[S1]+E[S2])/1000 < 1.
S1_mean = exp(1 + 0.5);           % E[S1] = exp(mu + sigma^2/2) with mu=sigma=1
S2_mean = 0.75;                   % E[S2] = (0.5+1)/2
lam_stab = 1000 / (S1_mean + S2_mean);
printf('\nA.2 stability bound: lambda < %.2f req/sec\n', lam_stab);

lambdas_A2 = [50];                % put all five values here: [50 100 150 200 250]
runs_A2    = cell(1, length(lambdas_A2));
for k = 1:length(lambdas_A2)
    printf('  running lambda = %d ... ', lambdas_A2(k)); fflush(stdout);
    runs_A2{k} = sim_webserver(lambdas_A2(k), 1e4);
    printf('done\n');
end
% Each runs_A2{k} has the usual fields (log_t, log_q1, log_q2, ...), so
% you can plot with e.g.:
%   plot(runs_A2{k}.log_t, runs_A2{k}.log_q1 + runs_A2{k}.log_q2)

% A.3: 95% CIs from 30 independent replications.
% Each run gives one time-weighted average of the number of type-i jobs.
% Those 30 values feed calc_ci (t-distribution for the mean, order
% statistics for the median). The script prints full-trace CIs below.
%
% What you need to do: plot MQS(t) for one run at each lambda, pick a
% cutoff, then redo the CIs after discarding everything before the
% cutoff. Compare with the full-trace CIs.

n_runs    = 30;
n_req_rep = 1e4;
lambdas   = [100, 180];

CI1_mean_all = zeros(length(lambdas), 2);
CI1_med_all  = zeros(length(lambdas), 2);
CI2_mean_all = zeros(length(lambdas), 2);
CI2_med_all  = zeros(length(lambdas), 2);

for li = 1:length(lambdas)
    lam = lambdas(li);
    N1 = zeros(n_runs, 1);
    N2 = zeros(n_runs, 1);
    % N1_tr = zeros(n_runs, 1);    % for your transient-removal CIs
    % N2_tr = zeros(n_runs, 1);

    printf('\nlambda = %d, running %d replications ', lam, n_runs);
    for k = 1:n_runs
        rk = sim_webserver(lam, n_req_rep);

        % full trace (includes the warmup)
        N1(k) = tw_mean(rk.log_q1, rk.log_t);
        N2(k) = tw_mean(rk.log_q2, rk.log_t);

        % pick a cutoff from the MQS(t) plot and compute the per-run
        % averages after removing the transient, e.g.:
        %   start = floor(0.3 * length(rk.log_t)) + 1;
        %   N1_tr(k) = tw_mean(rk.log_q1(start:end), rk.log_t(start:end));
        %   N2_tr(k) = tw_mean(rk.log_q2(start:end), rk.log_t(start:end));

        printf('.'); fflush(stdout);
    end
    printf(' done\n');

    [ci1_mean, ci1_med] = calc_ci(N1);
    [ci2_mean, ci2_med] = calc_ci(N2);
    CI1_mean_all(li,:) = ci1_mean;  CI1_med_all(li,:) = ci1_med;
    CI2_mean_all(li,:) = ci2_mean;  CI2_med_all(li,:) = ci2_med;

    printf('  full trace, type 1: mean CI [%.4f, %.4f], median CI [%.4f, %.4f]\n', ...
           ci1_mean(1), ci1_mean(2), ci1_med(1), ci1_med(2));
    printf('  full trace, type 2: mean CI [%.4f, %.4f], median CI [%.4f, %.4f]\n', ...
           ci2_mean(1), ci2_mean(2), ci2_med(1), ci2_med(2));

    % Once you have N1_tr and N2_tr:
    %   [ci1m, ci1md] = calc_ci(N1_tr);
    %   [ci2m, ci2md] = calc_ci(N2_tr);
end

% A.4: Little's Law.
% Verify that lambda * R is close to N at lambda = 100, with and without
% transient removal. You need:
%
%   R = mean(r.resp_t1) + mean(r.resp_t2)      end-to-end response time (ms)
%   N = tw_mean(r.log_q1 + r.log_q2, r.log_t)  time-weighted total jobs
%
% Watch the units. lambda is in req/sec and R is in ms, so convert lambda
% to req/ms before comparing:
%
%   lambda_ms = lambda / 1000;
%   N_little  = lambda_ms * R;      % should be close to N
%
% Then repeat after removing the transient and compare the two errors.

% Summary of everything the script printed.
% If the command window has scrolled, the numbers are reprinted below.

printf('\nSummary\n');

printf('\nA.1  (lambda = %d, %d requests, t_end = %.1f ms)\n', lambda, max_req, r.t_end);
printf('  R1 = %.4f ms,  R2 = %.4f ms,  R end-to-end = %.4f ms\n', ...
       mean(r.resp_t1), mean(r.resp_t2), mean(r.resp_t1)+mean(r.resp_t2));
printf('  jobs served per second: type-1 = %.4f, type-2 = %.4f\n', tp1, tp2);
printf('  mean in system (tw):  type-1 = %.4f,  type-2 = %.4f\n', N1_single, N2_single);

printf('\nA.2  stability bound: lambda < %.2f req/sec\n', lam_stab);

printf('\nA.3  full-trace CIs from %d replications, 95%% level:\n', n_runs);
for li = 1:length(lambdas)
    printf('  lambda = %d\n', lambdas(li));
    printf('    type 1: mean CI [%.4f, %.4f], median CI [%.4f, %.4f]\n', ...
           CI1_mean_all(li,1), CI1_mean_all(li,2), CI1_med_all(li,1), CI1_med_all(li,2));
    printf('    type 2: mean CI [%.4f, %.4f], median CI [%.4f, %.4f]\n', ...
           CI2_mean_all(li,1), CI2_mean_all(li,2), CI2_med_all(li,1), CI2_med_all(li,2));
end
