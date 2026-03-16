pkg load statistics
clc;
clear all;
close all;

# ============================================================
# Lab 1 - Part 2
# Exponential Distribution
# ============================================================

# We use the vector suggested in the assignment
k = 0:0.00001:8;

# The given mean values are 1/lambda = [0.5, 1.5, 4]
mean_values = [0.5, 1.5, 4];
lambda = 1 ./ mean_values;

# ------------------------------------------------------------
# TASK A
# PDF of three exponential distributions
# exppdf(x, mu) in Octave takes the mean value mu = 1/lambda
# ------------------------------------------------------------
pdf_exp = zeros(columns(mean_values), columns(k));

for i = 1:columns(mean_values)
  pdf_exp(i, :) = exppdf(k, mean_values(i));
endfor

colors = "rgb";

figure(1);
hold on;

for i = 1:columns(mean_values)
  plot(k, pdf_exp(i, :), colors(i), "linewidth", 1.5);
endfor

hold off;
grid on;
title("Probability Density Function of Exponential Distribution");
xlabel("k values");
ylabel("Probability values");
legend("mean = 0.5 sec (lambda = 2)", ...
       "mean = 1.5 sec (lambda = 2/3)", ...
       "mean = 4 sec (lambda = 1/4)");


# ------------------------------------------------------------
# TASK B
# CDFs of X1, X2, X3 and X = min(X1, X2, X3)
# If Xi ~ Exp(lambda_i), then min(X1, X2, X3) ~ Exp(lambda1 + lambda2 + lambda3)
# ------------------------------------------------------------
cdf_exp = zeros(columns(mean_values), columns(k));

for i = 1:columns(mean_values)
  cdf_exp(i, :) = expcdf(k, mean_values(i));
endfor

lambda_sum = sum(lambda);
mean_min = 1 / lambda_sum;

cdf_min = expcdf(k, mean_min);

figure(2);
hold on;

for i = 1:columns(mean_values)
  plot(k, cdf_exp(i, :), colors(i), "linewidth", 1.5);
endfor

plot(k, cdf_min, "k", "linewidth", 1.8);

hold off;
grid on;
title("Exponential Cumulative Distribution Function for X1, X2, X3 and min(X1,X2,X3)");
xlabel("k values");
ylabel("Probability values");
legend("X1: mean = 0.5 sec", ...
       "X2: mean = 1.5 sec", ...
       "X3: mean = 4 sec", ...
       "X = min(X1,X2,X3)");

fprintf("Rate of X = min(X1, X2, X3): lambda = %.6f\n", lambda_sum);
fprintf("Mean value of X = min(X1, X2, X3): %.6f sec\n\n", mean_min);


# ------------------------------------------------------------
# TASK C
# Memoryless property of the exponential distribution
# X ~ Exp with mean 1/lambda = 2.5 sec
# ------------------------------------------------------------
mean_X = 2.5;

index_20000 = 20000;
index_25000 = 25000;
index_45000 = 45000;

t_20000 = k(index_20000);
t_25000 = k(index_25000);
t_45000 = k(index_45000);

Prob_1 = 1 - expcdf(t_20000, mean_X);
Prob_2 = (1 - expcdf(t_45000, mean_X)) / (1 - expcdf(t_25000, mean_X));

# Direct memoryless check: P(X > t + s | X > s) = P(X > t)
Prob_memoryless_check = 1 - expcdf(t_45000 - t_25000, mean_X);

fprintf("For X ~ Exp with mean = 2.5 sec:\n\n");

fprintf("Index 20000 corresponds to time t = %.5f sec\n", t_20000);
fprintf("Index 25000 corresponds to time s = %.5f sec\n", t_25000);
fprintf("Index 45000 corresponds to time t+s = %.5f sec\n\n", t_45000);

fprintf("Pr(X > 20000) = %.6f\n", Prob_1);
fprintf("Pr(X > 45000 | X > 25000) = %.6f\n", Prob_2);
fprintf("Memoryless check Pr(X > (45000-25000)) = %.6f\n", Prob_memoryless_check);
