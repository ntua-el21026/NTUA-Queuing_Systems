pkg load statistics
clc;
clear all;
close all;

# TASK_A: PMF of Poisson distributions

# Define the lambda parameters
lambda = [7,20,40,70];
# Define the k values
k = 0:100;

# Define the Poisson PMF values
for i=1:columns(lambda)
  poisson(i,:) = poisspdf(k,lambda(i));
endfor

colors = "rbkm";
figure(1);
hold on;

# Plotting the Poisson PMF
for i=1:columns(lambda)
  stem(k,poisson(i,:),colors(i),"linewidth",1.2);
endfor
hold off;
title("Probability Mass Function of the Poisson distribution");
xlabel("k values");
ylabel("Probability values");
# legend for each case
legend1 = sprintf('lambda = %d',lambda(1));
legend2 = sprintf('lambda = %d',lambda(2));
legend3 = sprintf('lambda = %d',lambda(3));
legend4 = sprintf('lambda = %d',lambda(4));
legend(legend1,legend2,legend3,legend4)


# TASK_B: Binomial approximation to Poisson

# We want Poisson(40)
lambda_constant = 40;

# Define the n and p parameters for the approximation
n = [40,120,360,1080,40000];
p = lambda_constant./n;

# Define the Binomial PMF values
for i=1:columns(n)
  binomial(i,:) = binopdf(k,n(i),p(i));
endfor

colors = "ygmbr";
figure(2);
subplot(2,1,1);
hold on;

# Plotting the Bino(n,p) for all the values of n
for i=1:columns(n)
  stem(k,binomial(i,:),colors(i),"linewidth",1.2);
endfor
title("PMF of the Binomial distribution");
xlabel("k values");
ylabel("Probability values");
legend("n = 40  &  p = 1","n = 120  &  p = 1/3 ≈ 0.333","n = 360  &  p = 1/9 ≈ 0.111","n = 1080  &  p = 1/27 ≈ 0.037","n = 40000  &  p = 1/1000 = 0.001");
hold off;

subplot(2,1,2);
hold on;

# In order to zoom in and notice the approximation
# better, we will work in the following way:

# Obtain the position of the Po(40) from the 1st Fig.
index = find(lambda == 40);

# Plot the Bino(n,p) without the first (n=40) value
for i=2:columns(n)
  stem(k,binomial(i,:),colors(i),"linewidth",1.2);
endfor

# Include the Po(40) from the 1st Fig. with the same color
stem(k,poisson(index,:),"k","linewidth",1.2);
hold off;
title("Poisson as the limit of the Binomial distribution");
xlabel("k values");
ylabel("Probability values");
legend("n = 120  &  p = 1/3 ≈ 0.333","n = 360  &  p = 1/9 ≈ 0.111","n = 1080  &  p = 1/27 ≈ 0.037","n = 40000  &  p = 1/1000 = 0.001","Po(40)");


# TASK_C:
# Regarding the Po(40) distribution,
# compute the mean value and variance

# for the mean (or the "expected") value, we will use the definition
mean_value_pois = 0;
for i=0:(columns(poisson(index,:))-1)
  mean_value_pois = mean_value_pois + i.*poisson(index,i+1);
endfor
display("Mean value of Poisson with parameter lambda = 40 is");
display(mean_value_pois);

# for the variance, we will use the moments formula
second_moment_pois = 0;
for i=0:(columns(poisson(index,:))-1)
  second_moment_pois = second_moment_pois + i.*i.*poisson(index,i+1);
endfor
variance_pois = second_moment_pois - mean_value_pois.^2;
display("Variance of Poisson with parameter lambda = 40 is");
display(variance_pois);

# Mean value of the Binomial
mean_value_bino = zeros(1,columns(n));
for index_b=1:columns(n)
  for i=0:(columns(binomial(index_b,:))-1)
    mean_value_bino(1,index_b) = mean_value_bino(1,index_b) + i.*binomial(index_b,i+1);
  endfor
endfor
display("Mean value of the Binomial with parameter n and p");
display(mean_value_bino);

# Variance of the Binomial
variance_bino = zeros(1,columns(n));
second_moment_bino = zeros(1,columns(n));
for index_b=1:columns(n)
  for i=0:(columns(binomial(index_b,:))-1)
    second_moment_bino(1,index_b) = second_moment_bino(1,index_b) + i.*i.*binomial(index_b,i+1);
  endfor
  variance_bino(1,index_b) = second_moment_bino(1,index_b) - (mean_value_bino(1,index_b))^2;
endfor
display("Variance of the Binomial as n increases");
display(variance_bino);


# TASK_D: Poisson decomposition/superposition.

# a) Plot in a common diagram the 3 distributions

# Define the lambda parameter
lambdaX = 30;

# Define the Poisson PMF values regarding the number of calls
PMF_X = poisspdf(k,lambdaX);

# Define the probabilities regarding the calls (external, internal)
p_ext = 0.3;
q_int = 1-p_ext;

# Thus, each type of call-event will follow a Poisson distribution

# Define the lambda parameters:
lambdaX_1 = lambdaX*p_ext;
lambdaX_2 = lambdaX*q_int;

# Define the PMF values for each type of call
PMF_X1 = poisspdf(k,lambdaX_1);
PMF_X2 = poisspdf(k,lambdaX_2);

colors = "rbg";
figure(3);
hold on;

# Plotting the Poisson calls
stem(k,PMF_X1,colors(1),"linewidth",1.2);
stem(k,PMF_X2,colors(2),"linewidth",1.2);
stem(k,PMF_X,colors(3),"linewidth",1.2);
hold off;
title("Probability Mass Function of the Poisson calls");
xlabel("k values");
ylabel("Probability values");
# legend for each case
legend1 = sprintf('lambdaX_1 = %d', lambdaX_1);
legend2 = sprintf('lambdaX_2 = %d', lambdaX_2);
legend3 = sprintf('lambdaX = %d', lambdaX);
legend(legend1,legend2,legend3)

# Calculation of probabilities

# using the number of combinations a.k.a. "binomial coefficient"
# (or directly using the binopdf from the statistics pkg)

# the probability (given n and k) can be found in either of the following ways:

# b)
# Define the number of n trials (number of calls made)
number_of_calls_b = 2;
# Define the k "successes", that, in our case, is
# the number of external calls with probability p_ext
k_ext_b = 2;
Prob_b = (nchoosek(number_of_calls_b,k_ext_b))*(p_ext^k_ext_b)*(q_int^(number_of_calls_b-k_ext_b));
# or directly with the use of binopdf (requires the statistics pkg)
Prob_b_bino = binopdf(k_ext_b,number_of_calls_b,p_ext);
display("The probability that both of the next two calls are external is");
display(Prob_b);
display("Check with binopdf");
display(Prob_b_bino);

# c)
number_of_calls_c = 5;
k_ext_c = 2;
Prob_c = (nchoosek(number_of_calls_c,k_ext_c))*(p_ext^k_ext_c)*(q_int^(number_of_calls_c-k_ext_c));
# or, with the use of binopdf,
Prob_c_bino = binopdf(k_ext_c,number_of_calls_c,p_ext);
display("The probability that out of the 5 calls, two calls exactly are external, is");
display(Prob_c);
display("Check with binopdf");
display(Prob_c_bino);
